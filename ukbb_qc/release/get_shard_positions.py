import logging
import os

import hail as hl

from tgg_methods.batch.batch_utils import (
    init_arg_parser,
    init_job,
    localize_file,
    run_batch,
)

from ukbb_qc.resources.basics import release_vcf_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("get_vcf_start_stop_pos")
logger.setLevel(logging.INFO)


def parse_args():
    """Parse command line args."""
    parser = init_arg_parser(default_billing_project="gnomad-production")
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "-d",
        "--docker",
        help="Docker image to use",
        default="gcr.io/broad-mpg-gnomad/tgg-methods-vm:20210623",
    )
    parser.add_argument(
        "-c",
        "--cpu",
        help="Number of CPUs to use for batch job. Default is 0.25.",
        default=0.25,
        type=float,
    )
    parser.add_argument(
        "--output-dir",
        help="Where to write output files.",
        default="gs://broad-ukbb/broad.freeze_6/temp",
    )
    args = parser.parse_args()
    return args


def main():
    """
    Get start and end position of VCF shards.

    Required for return of results to UKBB.
    """

    args = parse_args()
    logger.info("Getting VCF shards...")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)
    shards = hl.hadoop_ls(f"{release_vcf_path(*tranche_data)}/*.bgz")

    logger.info("Preparing to start batch job...")
    with run_batch(args, batch_name="UKBB vcf start stops") as batch:
        for s in shards:
            # Each shard looks like this:
            # {'path': 'gs://broad-ukbb/broad.freeze_6/release/vcf/sharded_vcf/broad.freeze_6.bgz/part-09999.bgz', 'size_bytes': 1666945357,
            # 'size': '1.6G', 'is_dir': False, 'modification_time': 'Mon Aug 17 22:24:28 EDT 2020', 'owner': 'kchao'}
            shard = s["path"]
            j = init_job(batch=batch, name=s, image=args.image, cpu=args.cpu)
            j.command(
                f"""gcloud -q auth activate-service-account --key-file=/gsa-key/key.json"""
            )
            local_vcf_shard = localize_file(j, shard)
            shard_name = os.path.basename(shard)
            j.command(
                f"""bcftools query -f '%CHROM %POS %REF %ALT\n' {local_vcf_shard} > temp.txt"""
            )
            j.command(f"""echo {shard_name} >> {shard_name}.txt""")
            j.command(f"""head -1 temp.txt >> {shard_name}.txt""")
            j.command(f"""tail -1 temp.txt >> {shard_name}.txt""")
            j.command(f"""gsutil -m cp {shard_name}.txt {args.output_dir}""")


if __name__ == "__main__":
    main()
