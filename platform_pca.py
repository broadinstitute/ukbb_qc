from gnomad_hail import *
from gnomad_hail.utils.generic import filter_to_autosomes
from ukbb_qc.resources import *
from ukbb_qc.utils import *
import hail as hl
import argparse
import numpy as np
import hdbscan


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("platform_pca")
logger.setLevel(logging.INFO)


# Note: The following functions are from Laurent's myoseq sample qc should probably move to a common location
def compute_callrate_mt(mt: hl.MatrixTable, intervals: hl.Table, reference_genome: str) -> hl.MatrixTable:
    callrate_mt = filter_to_autosomes(mt, reference_genome=reference_genome)
    callrate_mt = callrate_mt.annotate_rows(interval=intervals[callrate_mt.locus].target)
    callrate_mt = callrate_mt.filter_rows(hl.is_defined(callrate_mt.interval) & (hl.len(callrate_mt.alleles) == 2))
    callrate_mt = callrate_mt.select_entries(GT=hl.or_missing(hl.is_defined(callrate_mt.GT), hl.struct()))
    callrate_mt = callrate_mt.group_rows_by(callrate_mt.interval).aggregate(callrate=hl.agg.fraction(hl.is_defined(callrate_mt.GT)))
    return callrate_mt


def run_platform_pca(callrate_mt: hl.MatrixTable) -> Tuple[List[float], hl.Table, hl.Table]:
    callrate_mt = callrate_mt.annotate_entries(callrate=hl.int(callrate_mt.callrate > 0.25))
    
    # Center until Hail's PCA does it for you
    callrate_mt = callrate_mt.annotate_rows(mean_callrate=hl.agg.mean(callrate_mt.callrate))
    callrate_mt = callrate_mt.annotate_entries(callrate=callrate_mt.callrate - callrate_mt.mean_callrate)
    eigenvalues, scores, loadings = hl.pca(callrate_mt.callrate, compute_loadings=True)
    logger.info('Eigenvalues: {}'.format(eigenvalues))
    return eigenvalues, scores, loadings


def assign_platform_pcs(platform_pc_table: hl.Table, pc_scores_ann: str = 'scores', hdbscan_min_cluster_size: int = 500, hdbscan_min_samples: int = None) -> hl.Table:
    """
    Function assumes that platform_pc_table contains columns named 'combined_sample', 'gross_platform' (known labels), 'callratePC<n>'
    :param Table platform_pc_table: Table containing samples and callrate PCs
    :param str pc_scores_ann: annotation of the scores in the platform_pc_table
    :param int hdbscan_min_cluster_size: Minimum cluster size parameter for HDBSCAN
    :param int hdbscan_min_samples: Minimum samples parameter for HDBSCAN. If not specified, hdbscan_min_cluster_size is used.
    :return: Table containing samples, callrate PCs, and imputed platform labels
    :rtype: Table
    """
    # Read and format data for clustering
    data = platform_pc_table.to_pandas()
    callrate_data = np.matrix(data[pc_scores_ann].tolist())
    logger.info('Assigning platforms to {} exome samples in MT...'.format(len(callrate_data)))

    # Cluster data
    clusterer = hdbscan.HDBSCAN(min_cluster_size=hdbscan_min_cluster_size, min_samples=hdbscan_min_samples)
    cluster_labels = clusterer.fit_predict(callrate_data)
    n_clusters = len(set(cluster_labels)) - (-1 in cluster_labels)  # NOTE: -1 is the label for noisy (un-classifiable) data points
    logger.info('Found {} unique platforms during platform imputation...'.format(n_clusters))

    data['qc_platform'] = cluster_labels
    ht = hl.Table.from_pandas(data, key=[*platform_pc_table.key])
    ht = ht.annotate(qc_platform=hl.int32(ht.qc_platform))
    return ht


def main(args):
    hl.init(log='/platform_pca.log')

    data_source = args.data_source
    freeze = args.freeze

    if args.compute_callrate_mt:
        logger.info('Preparing data for platform PCA...')
        intervals = hl.import_locus_intervals(ukbb_calling_intervals_path, reference_genome='GRCh38')
        mt = get_ukbb_data(args.data_source, args.freeze, split=False, adj=True)
        mt = remove_hard_filter_samples(data_source, freeze, mt) 
        callrate_mt = compute_callrate_mt(mt, intervals, reference_genome='GRCh38')
        callrate_mt.write(callrate_mt_path(data_source, freeze), args.overwrite)

    if args.run_platform_pca:
        logger.info("Running platform PCA...")
        callrate_mt = hl.read_matrix_table(callrate_mt_path(data_source, freeze))
        eigenvalues, scores_ht, loadings_ht = run_platform_pca(callrate_mt)
        scores_ht.write(platform_pca_scores_ht_path(data_source, freeze), overwrite=args.overwrite)
        loadings_ht.write(platform_pca_loadings_ht_path(data_source, freeze), overwrite=args.overwrite)
        # Regeneron freeze 4 Eigenvalues: [26489244.935849957, 2039950.6985898241, 1407875.3058482022, 1082106.1507608977, 373810.0800184624, 361301.2291929654, 324435.7483132424, 205912.4810229146, 196196.71017912056, 159808.25367132248]

    if args.assign_platforms:
        logger.info("Assigning platforms based on platform PCA clustering")
        scores_ht = hl.read_table(platform_pca_scores_ht_path(data_source, freeze))
        platform_ht = assign_platform_pcs(scores_ht, hdbscan_min_cluster_size=args.hdbscan_min_cluster_size, hdbscan_min_samples=args.hdbscan_min_samples)
        platform_ht.write(platform_pca_results_ht_path(data_source, freeze), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)

    parser.add_argument('--compute_callrate_mt', help='Computes an interval by sample mt of callrate that will be used to compute platform PCs', action='store_true')
    parser.add_argument('--run_platform_pca', help='Runs platform PCA (assumes callrate MT was computed)', action='store_true')
    parser.add_argument('--assign_platforms', help='Assigns platforms based on callrate PCA results using HDBSCAN', action='store_true')
    parser.add_argument('--hdbscan_min_samples', help='Minimum samples parameter for HDBSCAN. If not specified, --hdbscan_min_cluster_size is used.', type=int, required=False)
    parser.add_argument('--hdbscan_min_cluster_size', help='Minimum cluster size parameter for HDBSCAN.', type=int, default=100)

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
