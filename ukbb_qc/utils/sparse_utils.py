import logging
from gnomad_hail.utils.generic import bi_allelic_expr, filter_to_autosomes
from ukbb_qc.resources.sample_qc import interval_qc_path


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sparse_utils")
logger.setLevel(logging.INFO)


def compute_callrate_dp_mt(
    data_source: str,
    freeze: int,
    mt: hl.MatrixTable,
    intervals_ht: hl.Table,
    bi_allelic_only: bool = True,
    autosomes_only: bool = True,
    match: bool = True,
) -> None:
    """
    Computes sample n_defined and n_count (inputs to callrate) and mean depth per interval.
    Mean depth annotation is used during interval QC.
    Callrate annotations are used during platform PCA.
    Writes call rate mt (Grouped MatrixTable) keyed by intervals row-wise and samples column-wise.
    NOTE: This function requires a densify! Please use an autoscaling cluster.

    :param data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param MatrixTable mt: Input MatrixTable.
    :param Table intervals_ht: Table with capture intervals relevant to input MatrixTable.
    :param bi_allelic_only: If set, only bi-allelic sites are used for the computation.
    :param autosomes_only: If set, only autosomal intervals are used.
    :param matches: If set, returns all intervals in interval_ht that overlap the locus in the input MT.
    :return: None
    """
    logger.info("Computing call rate and mean DP MatrixTable")
    if len(intervals_ht.key) != 1 or not isinstance(
        intervals_ht.key[0], hl.expr.IntervalExpression
    ):
        logger.warning(
            f"Call rate matrix computation expects `intervals_ht` with a key of type Interval. Found: {intervals_ht.key}"
        )
    logger.info("Densifying...")
    mt = hl.experimental.densify(mt)

    logger.info(
        "Filtering out lines that are only reference or not covered in capture intervals"
    )
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    intervals_ht = intervals_ht.annotate_rows(interval_label=intervals_ht.interval)
    mt = mt.annotate_rows(
        interval=intervals_ht.index(mt.locus, all_matches=match).interval_label
    )
    mt = mt.filter_rows(hl.is_defined(mt.interval))

    if autosomes_only:
        mt = filter_to_autosomes(mt)

    if bi_allelic_only:
        mt = mt.filter_rows(bi_allelic_expr(mt))

    if match:
        mt = mt.explode_rows("interval")

    logger.info(
        "Grouping MT by interval and calculating n_defined, n_count, and mean_dp..."
    )
    mt = mt.group_rows_by(**_mt._interval_key).aggregate(
        n_defined=hl.agg.count_where(hl.is_defined(mt.GT)),
        count=hl.agg.count(),
        mean_dp=hl.agg.mean(hl.if_else(hl.is_defined(mt.DP), mt.DP, 0)),
    )
    mt.write(callrate_mt_path(data_source, freeze, interval_filtered=False))


def get_sparse_sample_qc(mt: hl.MatrixTable, autosomes_only: bool = True) -> hl.Table:
    """
    Calculates approximation for call rate and mean DP on a sparse MatrixTable. 

    :param MatrixTable mt: Input sparse MatrixTable
    :param bool autosomes_only: Whether to use interval QC results calculated on autosomes only.
    :return: Table of samples with their call rate and mean DP
    :rtype: hl.Table
    """
    # Filter MatrixTable to autosomes if set
    if autosomes_only:
        mt = filter_to_autosomes(mt)
        chrom = "autosomes"
    else:
        chrom = None

    # Read in interval QC table and get total number of bases in high coverage intervals
    interval_qc_ht = hl.read_table(interval_qc_path(data_source, freeze, chrom=chrom))
    good_intervals_ht = interval_qc_ht.filter(
        interval_qc_ht[cov_filter_field] > pct_samples
    )
    good_intervals_ht = good_intervals_ht.annotate(
        interval_label=good_intervals_ht.interval
    )
    good_intervals = good_intervals_ht.aggregate(
        hl.agg.collect(good_intervals_ht.interval)
    )
    total_bases = good_intervals_ht.aggregate(
        hl.agg.sum(
            hl.abs(
                good_intervals_ht.interval.end.position
                - good_intervals_ht.interval.start.position
            )
        )
    )
    logger.info(f"Total number of bases in high coverage intervals: {total_bases}")

    # Add END entry to rows
    mt = mt.annotate_rows(
        end=hl.or_missing(hl.len(mt.alleles) == 1, hl.agg.max(mt.END))
    )

    # Filter MatrixTable to high coverage intervals
    logger.info(f"Total number of rows in input MatrixTable: {mt.count_rows()}")
    mt = mt.annotate_rows(
        interval_qc_pass=(
            hl.is_defined(good_intervals_ht[mt.row_key])
            | hl.is_defined(good_intervals_ht[hl.locus(mt.locus.contig, mt.end)])
        )
    )
    mt = mt.annotate_rows(interval=good_intervals_ht[mt.row_key].interval_label)
    logger.info(
        f"Number of rows after filtering to high coverage intervals: {mt.count_rows()}"
    )

    # Add reference block start annotation to rows
    # This is to only count bases in reference blocks that overlap high coverage intervals
    mt = mt.annotate_rows(
        start=hl.case()
        .when(
            hl.len(mt.alleles) == 1,
            hl.if_else(
                mt.locus.position > mt.interval.start.position,
                mt.locus.position,
                mt.interval.start.position,
            ),
        )
        .default(mt.locus.position)
    )

    # Add reference block interval annotation to entries
    # Uses row annotation start for start of interval
    # Adjusts end of interval to the minimum value between END entry and end of interval
    mt = mt.annotate_entries(
        ref_block_interval=(
            hl.or_missing(
                mt.LGT.is_hom_ref(),
                hl.parse_locus_interval(
                    hl.format(
                        "%s:%s-%s",
                        mt.locus.contig,
                        mt.start,
                        hl.if_else(
                            mt.END < mt.interval.end.position,
                            mt.END,
                            mt.interval.end.position,
                        ),
                    )
                ),
            )
        )
    )

    # Add reference block length annotation to entries
    mt = mt.annotate_entries(
        ref_block_len=(
            hl.or_missing(
                mt.LGT.is_hom_ref(),
                mt.ref_block_interval.end.position
                - mt.ref_block_interval.start.position,
            )
        )
    )

    # Use ref block length entry to calculate mean DP
    mt = mt.annotate_cols(
        mean_dp=hl.agg.sum(
            hl.if_else(mt.LGT.is_hom_ref(), mt.DP * mt.ref_block_len, mt.DP)
        )
        / total_bases
    )

    # Use ref block length entry to create n_called column annotation
    # n_called is the total number of bases in reference blocks and non ref sites for a sample
    mt = mt.annotate_cols(
        n_nonref=(hl.agg.count_where(mt.LGT.is_non_ref())),
        total_ref_block_len=hl.agg.sum(mt.ref_block_len),
    )
    ht = mt.annotate_cols(n_called=mt.n_nonref + mt.total_ref_block_len).cols()

    # Calculate callrate and return
    ht = ht.annotate(callrate=(ht.total_ref_block_len + ht.n_nonref) / total_bases)
    return ht
