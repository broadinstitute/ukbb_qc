from ukbb_qc.resources.sample_qc import interval_qc_path


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
