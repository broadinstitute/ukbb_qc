from gnomad_hail import *
from ukbb_qc.resources import *


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("rank_rf")
logger.setLevel(logging.INFO)


def run_sanity_checks(ht: hl.Table) -> None:
    """
    Runs and prints sanity checks on rank table.

    :param Table ht: input ranks Table
    :return: Nothing
    :rtype: None
    """
    print(ht.aggregate(hl.struct(was_split=hl.agg.counter(ht.was_split),
                                 has_biallelic_rank=hl.agg.counter(hl.is_defined(ht.biallelic_rank)),
                                 was_singleton=hl.agg.counter(ht.singleton),
                                 has_singleton_rank=hl.agg.counter(hl.is_defined(ht.singleton_rank)),
                                 was_split_singleton=hl.agg.counter(ht.singleton & ~ht.was_split),
                                 has_biallelic_singleton_rank=hl.agg.counter(hl.is_defined(ht.biallelic_singleton_rank)))))


def create_rf_rank(data_source: str, freeze: int, run_hash: str) -> None:
    """
    Creates a ranked table for a RF run and writes it to its correct location in annotations.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str run_hash: RF run hash
    :return: Nothing
    :rtype: None
    """
    logger.info(f"Creating rank file for {data_source}.freeze_{freeze} RF run {run_hash}")

    if not hl.hadoop_exists(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/rf_{run_hash}.ht/_SUCCESS'):
        adj_ht = add_adj_annotations(data_source, freeze)
        ht = hl.read_table(rf_path(data_source, freeze, 'rf_result', run_hash=run_hash))
        ht = ht.annotate(**adj_ht[ht.key],
                         score=ht.rf_probability['TP'])

        # Write to temp location as result will be overwritten
        ht.write(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/rf_{run_hash}.ht', overwrite=True)
    ht = hl.read_table(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/rf_{run_hash}.ht')

    ht = add_rank(ht,
                  score_expr=1-ht.score,
                  subrank_expr={
                      'singleton_rank': ht.singleton,
                      'biallelic_rank': ~ht.was_split,
                      'biallelic_singleton_rank': ~ht.was_split & ht.singleton,
                      'adj_rank': ht.ac_adj > 0,
                      'adj_biallelic_rank': ~ht.was_split & (ht.ac_adj > 0),
                      'adj_singleton_rank': ht.singleton & (ht.ac_adj > 0),
                      'adj_biallelic_singleton_rank': ~ht.was_split & ht.singleton & (ht.ac_adj > 0)
                  }
                  )
    ht.write(rf_path(data_source, freeze, 'rf_result', run_hash=run_hash), overwrite=True)


def create_vqsr_rank_ht(data_source: str, freeze: int):
    """
    Creates a rank table for VQSR and writes it to its correct location in annotations.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Nothing
    :rtype: None
    """
    logger.info(f"Creating rank file for {data_source}.freeze_{freeze} VQSR")
    if not hl.utils.hadoop_exists(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/vqsr.ht/_SUCCESS'):
        rf_ht = hl.read_table(rf_annotated_path(data_source, freeze))
        adj_ht = add_adj_annotations(data_source, freeze)
        ht = hl.read_table(var_annotations_ht_path(data_source, freeze, 'vqsr')).repartition(1000)
        logger.info('Filtering to high_quality samples and n_nonref==1...')
        ht = ht.annotate(**adj_ht[ht.key],
                         **rf_ht[ht.key],
                         score=ht.info.VQSLOD,
                         culprit=ht.info.culprit)
        ht = ht.filter(ht.n_nonref > 0)
        # ht = ht.repartition(1000, shuffle=False)
        ht.write(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/vqsr.ht', overwrite=True)
    ht = hl.read_table(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/vqsr.ht')

    ht = add_rank(ht,
                  score_expr=-1*ht.score,
                  subrank_expr={
                      'singleton_rank': ht.singleton,
                      'biallelic_rank': ~ht.was_split,
                      'biallelic_singleton_rank': ~ht.was_split & ht.singleton,
                      'adj_rank': ht.ac_adj > 0,
                      'adj_biallelic_rank': ~ht.was_split & (ht.ac_adj > 0),
                      'adj_singleton_rank': ht.singleton & (ht.ac_adj > 0),
                      'adj_biallelic_singleton_rank': ~ht.was_split & ht.singleton & (ht.ac_adj > 0)
                  }
                  )
    ht.write(score_ranking_path(data_source, freeze, 'vqsr'), overwrite=True)


def add_adj_annotations(data_source: str, freeze: int) -> hl.Table:
    """
    Gets a Table with gnomad annotations needed to create ranked tables.
    Because this is computationally expensive, the table is written the first time it is created and then loaded after that.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Table with annotations
    :rtype: Table
    """
    if not hl.utils.hadoop_exists(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/adj_annotations.ht'):
        mt = get_ukbb_data(data_source, freeze, adj=True)
        ht = mt.annotate_rows(ac_adj=hl.agg.sum(mt.GT.n_alt_alleles())).rows().select('ac_adj')
        ht.repartition(1000, shuffle=False).write(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/adj_annotations.ht')
    return hl.read_table(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/adj_annotations.ht')


def create_binned_data(ht: hl.Table, data: str, data_source: str, freeze: int, n_bins: int) -> hl.Table:
    """
    Creates binned data from a rank Table grouped by rank_id (rank, biallelic, etc.), contig, snv, bi_allelic and singleton
    containing the information needed for evaluation plots.

    :param Table ht: Input rank table
    :param str data: Which data/run hash is being created
    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param int n_bins: Number of bins.
    :return: Binned Table
    :rtype: Table
    """

    # Count variants for ranking
    count_expr = {x: hl.agg.filter(hl.is_defined(ht[x]), hl.agg.counter(hl.cond(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv', 'indel'))) for x in ht.row if x.endswith('rank')}
    rank_variant_counts = ht.aggregate(hl.Struct(**count_expr))
    logger.info(f"Found the following variant counts:\n {pformat(rank_variant_counts)}")
    ht = ht.annotate_globals(rank_variant_counts=rank_variant_counts)

    # Load external evaluation data
    # clinvar_ht = hl.read_table(clinvar_ht_path)  # TODO: swap out for official table
    if not hl.utils.hadoop_exists('gs://broad-ukbb/resources/clinvar_20190609.pathogenic.ht'):
        clinvar_ht = hl.read_table('gs://gnomad/tmp/clinvar_20190609.ht')

        logger.info(f"Found {clinvar_ht.count()} variants in Clinvar before filtering")
        no_star_assertions = hl.literal({"no_assertion_provided", "no_assertion_criteria_provided", "no_interpretation_for_the_single_variant"})
        clinvar_ht = clinvar_ht.filter(hl.set(clinvar_ht.info.CLNREVSTAT).intersection(no_star_assertions).length() > 0, keep=False)
        logger.info(f"Found {clinvar_ht.count()} variants after removing variants without assertions")
        clinvar_ht = clinvar_ht.filter(clinvar_ht.info.CLNSIG.map(lambda x: x.lower()).map(lambda x: x.contains('pathogenic')).any(lambda x: x), keep=True)
        logger.info(f"Found {clinvar_ht.count()} variants after filtering to (likely) pathogenic variants")
        clinvar_ht = clinvar_ht.filter(hl.is_defined(clinvar_ht.info.CLNSIGCONF), keep=False)
        logger.info(f"Found {clinvar_ht.count()} variants after filtering to variants without CLNSIGCONF (conflicting clinical interpretations)")
        clinvar_ht.write('gs://broad-ukbb/resources/clinvar_20190609.pathogenic.ht')
    clinvar_ht = hl.read_table('gs://broad-ukbb/resources/clinvar_20190609.pathogenic.ht')

    # denovo_ht = get_validated_denovos_ht()
    # if data_type == 'exomes':
    #     denovo_ht = denovo_ht.filter(denovo_ht.gnomad_exomes.high_quality)
    # else:
    #     denovo_ht = denovo_ht.filter(denovo_ht.gnomad_genomes.high_quality)
    # denovo_ht = denovo_ht.select(validated_denovo=denovo_ht.validated, high_confidence_denovo=denovo_ht.Confidence == 'HIGH')
    ht_truth_data = hl.read_table(var_annotations_ht_path(data_source, freeze, 'truth_data'))
    fam_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, 'family_stats'))
    fam_ht = fam_ht.select(
        family_stats=fam_ht.family_stats[0]
    )
    ukbb_ht = get_ukbb_data(data_source, freeze).rows()
    ukbb_ht = ukbb_ht.select(
        fail_hard_filters=(ukbb_ht.info.QD < 2) | (ukbb_ht.info.FS > 60) | (ukbb_ht.info.MQ < 30)
    )

    ht = ht.annotate(
        **ht_truth_data[ht.key],
        **fam_ht[ht.key],
        **ukbb_ht[ht.key],
        # **denovo_ht[ht.key],
        clinvar=hl.is_defined(clinvar_ht[ht.key]),
        indel_length=hl.abs(ht.alleles[0].length()-ht.alleles[1].length()),
        rank_bins=hl.array(
            [hl.Struct(
                rank_id=rank_name,
                bin=hl.int(hl.ceil(hl.float(ht[rank_name] + 1) / hl.floor(ht.globals.rank_variant_counts[rank_name][hl.cond(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv', 'indel')] / n_bins)))
            )
                for rank_name in rank_variant_counts]
        )
    )

    ht = ht.explode(ht.rank_bins)
    ht = ht.transmute(
        rank_id=ht.rank_bins.rank_id,
        bin=ht.rank_bins.bin
    )
    ht = ht.filter(hl.is_defined(ht.bin))

    ht.write(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/score_binning.{data}.ht', overwrite=True)

    # Create binned data
    ht = hl.read_table(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/score_binning.{data}.ht')

    return (
        ht.group_by(
            rank_id=ht.rank_id,
            contig=ht.locus.contig,
            snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
            bi_allelic=hl.is_defined(ht.biallelic_rank),
            singleton=ht.singleton,
            release_adj=ht.ac_adj>0,
            bin=ht.bin
        )._set_buffer_size(20000)
        .aggregate(
            min_score=hl.agg.min(ht.score),
            max_score=hl.agg.max(ht.score),
            n=hl.agg.count(),
            n_ins=hl.agg.count_where(hl.is_insertion(ht.alleles[0], ht.alleles[1])),
            n_del=hl.agg.count_where(hl.is_deletion(ht.alleles[0], ht.alleles[1])),
            n_ti=hl.agg.count_where(hl.is_transition(ht.alleles[0], ht.alleles[1])),
            n_tv=hl.agg.count_where(hl.is_transversion(ht.alleles[0], ht.alleles[1])),
            n_1bp_indel=hl.agg.count_where(ht.indel_length == 1),
            n_mod3bp_indel=hl.agg.count_where((ht.indel_length % 3) == 0),
            n_clinvar=hl.agg.count_where(ht.clinvar),
            n_singleton=hl.agg.count_where(ht.singleton),
            # n_validated_de_novos=hl.agg.count_where(ht.validated_denovo),
            # n_high_confidence_de_novos=hl.agg.count_where(ht.high_confidence_denovo),
            n_de_novo=hl.agg.filter(ht.family_stats.unrelated_qc_callstats.AC[1] == 0, hl.agg.sum(ht.family_stats.mendel.errors)),
            n_trans_singletons=hl.agg.filter((ht.info_ac < 3) & (ht.family_stats.unrelated_qc_callstats.AC[1] == 1), hl.agg.sum(ht.family_stats.tdt.t)),
            n_untrans_singletons=hl.agg.filter((ht.info_ac < 3) & (ht.family_stats.unrelated_qc_callstats.AC[1] == 1), hl.agg.sum(ht.family_stats.tdt.u)),
            n_train_trans_singletons=hl.agg.count_where((ht.family_stats.unrelated_qc_callstats.AC[1] == 1) & (ht.family_stats.tdt.t == 1)),
            n_omni=hl.agg.count_where(ht.truth_data.omni),
            n_mills=hl.agg.count_where(ht.truth_data.mills),
            n_ukbb_array=hl.agg.count_where(ht.truth_data.ukbb_array),
            n_ukbb_array_con=hl.agg.count_where(ht.truth_data.ukbb_array_con),
            n_ukbb_array_con_common=hl.agg.count_where(ht.truth_data.ukbb_array_con_common),
            n_sib_singletons=hl.agg.count_where(ht.truth_data.sib_singletons),
            n_hapmap=hl.agg.count_where(ht.truth_data.hapmap),
            n_kgp_high_conf_snvs=hl.agg.count_where(ht.truth_data.kgp_high_conf_snvs),
            fail_hard_filters=hl.agg.count_where(ht.fail_hard_filters)
        )
    )


def main(args):

    def run_data(rank_func: Callable[..., None], rank_func_args: List[Any], score_type: str, data_source: str, freeze: int) -> None:
        """
        Wrapper for running script actions on given data.

        :param callable rank_func: Function for creating ranking file
        :param list of Any rank_func_args: Arguments to pass to the ranking function
        :param str score_type: Score being processed
        :param str data_source: 'regeneron' or 'broad'
        :param int freeze: One of the data freezes
        :return: Nothing -- this runs the script actions
        :rtype:  None
        """
        if score_type in ['vqsr']:
            rank_file_path = score_ranking_path(data_source, freeze, score_type)
        else:
            rank_file_path = rf_path(data_source, freeze, 'rf_result', run_hash=args.run_hash)

        if args.create_rank_file:
            rank_func(*rank_func_args)
        if args.run_sanity_checks:
            run_sanity_checks(hl.read_table(rank_file_path))
        if args.create_binned_file:
            ht = hl.read_table(rank_file_path)
            binned_ht = create_binned_data(ht, score_type, data_source, freeze, args.n_bins)
            binned_ht.write(score_ranking_path(data_source, freeze, score_type, binned=True), overwrite=True)

    hl.init(log='/create_rank.log')
    data_source = args.data_source
    freeze = args.freeze

    if args.run_hash:
        run_data(create_rf_rank, [data_source, freeze, args.run_hash], args.run_hash, data_source, freeze)

    if args.vqsr:
        run_data(create_vqsr_rank_ht, [data_source, freeze], 'vqsr', data_source, freeze)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)

    parser.add_argument('--run_hash', help='Run hash for RF results to be ranked.')
    parser.add_argument('--vqsr', help='When set, creates the VQSR rank file.', action='store_true')
    parser.add_argument('--create_rank_file', help='When set, creates ranking file.', action='store_true')
    parser.add_argument('--run_sanity_checks', help='When set, runs ranking sanity checks.', action='store_true')
    parser.add_argument('--create_binned_file', help='When set, creates binned ranked file.', action='store_true')
    parser.add_argument('--n_bins', help='Number of bins for the binned file (default: 100)', default=100, type=int)
    args = parser.parse_args()

    if not args.create_rank_file and not args.run_sanity_checks and not args.create_binned_file:
        sys.exit('Error: At least one of --create_rank_file, --run_sanity_checks or --create_binned_file must be specified.')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
