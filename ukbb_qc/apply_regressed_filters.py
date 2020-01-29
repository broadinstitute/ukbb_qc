import hail as hl


from gnomad_hail import *
from gnomad_hail.utils.sample_qc import *
from ukbb_qc.resources import * 


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('apply_regressed_filters')
logger.setLevel(logging.INFO)


def apply_regressed_filters(
        sample_qc_ht: hl.Table,
        pca_ht: hl.Table,
        overwrite: bool,
        filtering_qc_metrics: List[str],
        use_pc_square: bool = True,
        n_pcs: Optional[int] = None):

    sample_qc_ht = sample_qc_ht.annotate(**pca_ht[sample_qc_ht.s])
    sample_qc_ht.describe()

    sample_qc_ht = sample_qc_ht.repartition(100) # Repartition is needed for downstream operations
    sample_qc_ht = sample_qc_ht.checkpoint(
        'gs://broad-ukbb/broad.freeze_4/temp/qc.repartitioned.ht', 
        overwrite=overwrite, _read_if_exists=not overwrite)

    if n_pcs is None:
        n_pcs = sample_qc_ht.aggregate(hl.agg.min(hl.len(sample_qc_ht.scores)))

    logger.info("Computing regressed QC metrics filters using {} PCs for metrics: {}".format(
        n_pcs,
        ', '.join(filtering_qc_metrics)
    ))

    x_expr = [1.0] + [sample_qc_ht.scores[i] for i in range(0, n_pcs)]
    if use_pc_square:
        x_expr.extend([sample_qc_ht.scores[i]*sample_qc_ht.scores[i] for i in range(0, n_pcs)])

    lms = sample_qc_ht.aggregate(
        hl.struct(
            **{metric:
                hl.agg.linreg(
                    y=sample_qc_ht.sample_qc[metric],
                    x=x_expr
                 ) for metric in filtering_qc_metrics if not metric.startswith('log')
            }
        )
    )

    sample_qc_ht = sample_qc_ht.annotate_globals(
        lms=lms
    ).persist()

    def get_lm_prediction(metric: str):
        lm_pred_expr = sample_qc_ht.lms[metric].beta[0] + hl.sum(
                hl.range(n_pcs).map(lambda i: sample_qc_ht.lms[metric].beta[i + 1] * sample_qc_ht.scores[i])
        )
        if use_pc_square:
            lm_pred_expr = lm_pred_expr + hl.sum(
                hl.range(n_pcs).map(lambda i: sample_qc_ht.lms[metric].beta[i + n_pcs + 1] * sample_qc_ht.scores[i] * sample_qc_ht.scores[i])
        )
        return lm_pred_expr

    sample_qc_ht = sample_qc_ht.select(
        **{
            metric: sample_qc_ht.sample_qc[metric] - get_lm_prediction(metric)
            for metric in sample_qc_ht.lms
        }
    )

    sample_qc_ht = sample_qc_ht.persist()

    stratified_metrics_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        {metric: sample_qc_ht[metric] for metric in filtering_qc_metrics},
        metric_threshold={'n_singleton': (4.0, 8.0)}
    )

    sample_qc_ht = sample_qc_ht.annotate(
        **stratified_metrics_ht[sample_qc_ht.key]
    )
    sample_qc_ht = sample_qc_ht.annotate_globals(
        **stratified_metrics_ht.index_globals()
    )

    sample_qc_ht = sample_qc_ht.rename(
        {metric: f'{metric}_residual' for metric in filtering_qc_metrics}
    )

    sample_qc_ht.write('gs://broad-ukbb/broad.freeze_4/temp/regressed_metrics.ht', overwrite=overwrite)


data_source='broad'
freeze=4

apply_regressed_filters(
                        sample_qc_ht = hl.read_table(qc_temp_data_prefix(data_source, freeze) + 'outlier_sample_qc.ht'),
                        pca_ht = hl.read_table(ancestry_pca_scores_ht_path(data_source, freeze)),
                        overwrite=True, 
                        filtering_qc_metrics=['n_snp', 'r_ti_tv', 'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var'],
                        use_pc_square=True,
                        #use_pc_square=False, # NOTE: turning this to false filters a lot more samples
                        n_pcs=10
                        )
