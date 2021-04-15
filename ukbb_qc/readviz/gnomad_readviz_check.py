import argparse
import logging
import sqlite3
from typing import List


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("gnomad_readviz_check")
logger.setLevel(logging.INFO)


CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
"""
List of contigs to include.
"""

EXTRA_GENOME_CONTIGS = ["chrM"]
"""
Extra contigs to include for gnomAD genomes.

gnomAD exomes does not have chrM calls.
"""


def main(args):
    """
    This script extracts variants from gnomAD v2.1 (exomes) and v3.1 (genomes) that have fewer than the max number of samples available for readviz.

    This script outputs TSVs which will be imported into hail Tables in a downstream step.
    
    .. note::
        Script assumes that all sqlite databases are on the same machine that script is being run (i.e., local computer or cloud VM)
    """

    def _check_sample_counts(row_iter: sqlite3.Cursor) -> List[str]:
        """
        Iterate through sqlite curosr and extract variants with fewer than max number of samples displayed.

        :param sqlite3.Cursor row_iter: sqlite Cursor object containing variants to be checked.
        :return: List of variants that have fewer than the max number of samples displayed for readviz.
        :rtype: List[str]
        """
        lines = []

        for (
            chrom,
            pos,
            ref,
            alt,
            het_or_hom_or_hemi,
            n_available_samples,
        ) in row_iter:
            # Store variant if it has fewer than max number of samples displayed
            if n_available_samples < 3:
                lines.append(
                    "\t".join(
                        map(
                            str,
                            [
                                chrom,
                                pos,
                                ref,
                                alt,
                                het_or_hom_or_hemi,
                                n_available_samples,
                            ],
                        )
                    )
                    + "\n"
                )
        return lines

    if args.exomes:
        logger.info("Working on gnomAD v2.1 exomes...")
        # Fields in exomes databases are:
        # 'id', 'chrom', 'pos', 'ref', 'alt', 'het_or_hom_or_hemi',
        # 'variant_id', 'n_expected_samples', 'n_available_samples'
        # gnomAD v2.1 has 1000 sqlite databases per chromosome
        num_db_per_chrom = 1

        for contig in ["chr22"]:

            with open(f"{args.tsv_dir_path}/{contig}.tsv", "w") as o:
                logger.info("Writing header to file...")
                o.write(
                    "chrom\tpos\tref\talt\thet_or_hom_or_hemi\tn_available_samples\n"
                )

                # Check each database for contig
                for db in range(num_db_per_chrom):
                    # Add variable to format db name with correct number of zeroes
                    # e.g., correct name is combined_chr20_000.db, not combined_chr20_0.db
                    n_zeroes = 3 - len(str(db))

                    with sqlite3.connect(
                        f"{args.db_path}/combined_{contig}_{'0' * n_zeroes}{db}.db"
                    ) as conn:
                        row_iter = conn.execute(
                            "SELECT chrom, pos, ref, alt, het_or_hom_or_hemi, n_available_samples FROM t"
                        )
                        for line in _check_sample_counts(row_iter):
                            o.write(line)

    if args.genomes:
        logger.info("Working on gnomAD v3.1 genomes...")
        # Fields in genomes databases are:
        # 'id', 'chrom', 'pos', 'ref', 'alt', 'zygosity',
        # 'qual', 'combined_bamout_id', 'read_group_id'

        contigs = CONTIGS + EXTRA_GENOME_CONTIGS
        for contig in contigs:
            with open(f"{args.tsv_dir_path}/{contig}.tsv", "w") as o:
                logger.info("Writing header to file...")
                o.write(
                    "chrom\tpos\tref\talt\thet_or_hom_or_hemi\tn_available_samples\n"
                )

                # Base for gnomAD genomes databases is 'all_variants_s42811_gs50_gn857'
                with sqlite3.connect(
                    f"{args.db_path}/all_variants_s42811_gs50_gn857.{contig}.db"
                ) as conn:
                    # NOTE: Instead of storing 1 row per variant, het_or_hom_or_hemi, and the number of available tracks),
                    # the genome dbs store 1 row per available readviz track, so you have to group by to get the number of available samples
                    row_iter = conn.execute(
                        "SELECT chrom, pos, ref, alt, zygosity, count(*) FROM variants GROUP BY chrom, pos, ref, alt, zygosity"
                    )
                    for line in _check_sample_counts(row_iter):
                        o.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--exomes",
        help="Check for variants in gnomAD v2.1 exomes that do not have the max number of samples displayed (readviz). NOTE that gnomAD v2.1 is in GRCh37!",
        action="store_true",
    )
    parser.add_argument(
        "--genomes",
        help="Check for variants in gnomAD v3.1 genomes that do not have the max number of samples displayed (readviz)",
        action="store_true",
    )
    parser.add_argument(
        "--db_path",
        help="Absolute path to LOCAL directory that contains sqlite databases",
    )
    parser.add_argument(
        "--tsv_dir_path",
        help="Absolute path to LOCAL directory to store output (temporary) TSV files",
    )
    parser.add_argument(
        "--max_n_samples",
        type=int,
        help="Max number of samples to display per genotype",
        default=3,
    )

    args = parser.parse_args()
    main(args)
