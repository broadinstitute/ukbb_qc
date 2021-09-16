import logging
import os
import subprocess
from tqdm import tqdm

from tgg.batch.batch_utils import (
    localize_file,
    init_arg_parser,
    init_job,
    run_batch,
)

from ukbb_qc.resources.basics import release_vcf_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("repackage_vcf_shards")
logger.setLevel(logging.INFO)


DOCKER_IMAGE = "gcr.io/broad-mpg-gnomad/tgg-methods-vm:20210915"
"""
Docker image to be passed to hail batch. 

Must have bcftools, bgzip, tabix installed.
"""


def main(args):
    """
    Repackage each VCF shard so that the samples in the header are contained in a separate, 0-compressed bgzip block.

    Required for return of results to UKBB.

    .. note::
        Run this script in the same directory as `ukbb_header_reformat.sh`.
    """
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)
    output_dir = args.output_dir

    logger.info("Getting VCF shards...")
    path = f"{release_vcf_path(*tranche_data, contig=None)}/*bgz"
    shards = (
        subprocess.check_output(["gsutil", "ls", path])
        .decode("utf8")
        .strip()
        .split("\n")
    )

    # Process shards
    with run_batch(args, batch_name="Repackage VCF shards") as batch:
        for shard in tqdm(shards, unit="shards"):
            # Naming the output shard with 'repackaged' + gz extension since that is the extension DNANexus uses in their code
            output_name = f"{os.path.split(shard)[-1].split('.')[0]}.repackaged.gz"
            output_tabix = f"{output_name}.tbi"
            j = init_job(
                batch, f"{shard} repackage", DOCKER_IMAGE, args.cpu, args.memory,
            )
            j.command(
                f"""gcloud -q auth activate-service-account --key-file=/gsa-key/key.json"""
            )
            local_vcf_path = localize_file(j, shard, use_gcsfuse=True)
            ukbb_header_reformat_sh_contents = open("ukbb_header_reformat.sh").read()
            j.command(
                f"""cat <<"EOF" > ukbb_header_reformat.sh
{ukbb_header_reformat_sh_contents} 
EOF"""
            )
            j.command(
                f"""bash ukbb_header_reformat.sh {local_vcf_path} {output_name}"""
            )
            j.command(f"""gsutil -m cp {output_name} {output_dir}""")
            j.command(f"""gsutil -m cp {output_tabix} {output_dir}""")


if __name__ == "__main__":
    p = init_arg_parser(default_cpu=1, default_billing_project="gnomad-production")
    p.add_argument(
        "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    p.add_argument(
        "--output-dir",
        help="Where to write repackaged vcf output.",
        default="gs://broad-ukbb/broad.freeze_6/release/vcf/sharded_vcf/repackaged/",
    )
    args = p.parse_args()
    main(args)
