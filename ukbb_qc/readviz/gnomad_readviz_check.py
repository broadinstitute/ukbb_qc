import argparse
import logging
import sqlite3


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("gnomad_readviz_check")
logger.setLevel(logging.INFO)


CONTIGS = [i for i in range(23)] + ["X", "Y"]
"""
List of contigs to include.
"""

EXTRA_GENOME_CONTIGS = ["M"]
"""
Extra contigs to include for gnomAD genomes.

gnomAD exomes chrM calls.
"""


def main(args):
    """
    This script extracts variants from gnomAD v2.1 (exomes) and v3.1 (genomes) that have fewer than the max number of samples available for readviz.

    This script outputs TSVs which will be imported into hail Tables in a downstream step.
    
    Run this script locally.
    """
    if args.exomes:
        logger.info("Working on gnomAD v2.1 exomes...")
        # Fields in exomes databases are:
        # 'id', 'chrom', 'pos', 'ref', 'alt', 'het_or_hom_or_hemi',
        # 'variant_id', 'n_expected_samples', 'n_available_samples'
        # gnomAD v2.1 has 1000 sqlite databases per chromosome
        num_db_per_chrom = 1000

        variants = {}
        for contig in CONTIGS:
            variants[contig] = [
                "chrom\tpos\tref\talt\thet_or_hom_or_hemi\tn_available_samples\n"
            ]

            # Check each database for contig
            for db in range(num_db_per_chrom):
                # Add variable to format db name with correct number of zeroes
                # e.g., correct name is combined_chr20_000.db, not combined_chr20_0.db
                n_zeroes = 3 - len(str(db))

                with sqlite3.connect(
                    f"{args.db_path}/combined_chr{contig}_{'0' * n_zeroes}{db}.db"
                ) as conn:
                    row_iter = conn.execute(
                        "SELECT chrom, pos, ref, alt, het_or_hom_or_hemi, n_available_samples FROM t"
                    )
                    for (
                        chrom,
                        pos,
                        ref,
                        alt,
                        het_or_hom_or_hemi,
                        n_available_samples,
                    ) in row_iter:
                        # Store variant if it has fewer than max number of samples displayed
                        if n_available_samples < args.max_n_samples:
                            variants[contig].append(
                                f"{chrom}\t{pos}\t{ref}\t{alt}\t{het_or_hom_or_hemi}\t{n_available_samples}\n"
                            )

            # Write TSV with variants to output
            # NOTE: These TSVs will contain variants in GRCh37 --
            # they will need to be lifted over to GRCh38
            with open(f"{args.tsv_dir_path}/chr{contig}.tsv", "w") as o:
                for line in variants[contig]:
                    o.write(line)

    if args.genomes:
        logger.info("Working on gnomAD v3.1 genomes...")
        # Fields in genomes databases are:
        # 'id', 'chrom', 'pos', 'ref', 'alt', 'zygosity',
        # 'qual', 'combined_bamout_id', 'read_group_id'

        variants = {}
        contigs = CONTIGS + EXTRA_GENOME_CONTIGS
        for contig in contigs:
            variants[contig] = [
                "chrom\tpos\tref\talt\thet_or_hom_or_hemi\tn_available_samples\n"
            ]

            # Base for gnomAD genomes databases is 'all_variants_s42811_gs50_gn857'
            with sqlite3.connect(
                f"{args.db_path}/all_variants_s42811_gs50_gn857.chr{contig}.db"
            ) as conn:
                # NOTE: Instead of storing 1 row per variant, het_or_hom_or_hemi, and the number of available tracks),
                # the genome dbs store 1 row per available readviz track, so you have to group by to get the number of available samples
                row_iter = conn.execute(
                    "SELECT chrom, pos, ref, alt, zygosity, count(*) FROM variants GROUP BY chrom, pos, ref, alt, zygosity"
                )
                for (
                    chrom,
                    pos,
                    ref,
                    alt,
                    het_or_hom_or_hemi,
                    n_available_samples,
                ) in row_iter:
                    # Store variant if it has fewer than max number of samples displayed
                    if n_available_samples < args.max_n_samples:
                        variants[contig].append(
                            f"{chrom}\t{pos}\t{ref}\t{alt}\t{het_or_hom_or_hemi}\t{n_available_samples}\n"
                        )

            # Write TSV with variants to output
            with open(f"{args.tsv_dir_path}/chr{contig}.tsv", "w") as o:
                for line in variants[contig]:
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
