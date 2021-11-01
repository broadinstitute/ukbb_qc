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

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import parallel_file_exists

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
    coordinates_output_dir = "gs://broad-ukbb/broad.freeze_7/temp/coordinates"

    logger.info("Getting VCF shards...")
    contigs = []
    for i in range(1, 23):
        contigs.append(f"chr{i}")
    contigs.extend(["chrX", "chrY"])
    logger.info(f"Contigs: {contigs}")

    contig = args.contig
    if contig not in contigs:
        raise DataException(
            f"{contig} not in {contigs}. Please double check and restart!"
        )

    path = f"{release_vcf_path(*tranche_data, contig=contig)}/*bgz"
    shards = (
        subprocess.check_output(["gsutil", "ls", path])
        .decode("utf8")
        .strip()
        .split("\n")
    )
    logger.info("Found %i shards", len(shards))

    logger.info("Checking if any output files already exist...")
    output_coordinate_file_map = {}
    for shard in shards:
        output_name = f"ukb26041_{contig}_block{os.path.split(shard)[-1].split('.')[0].split('-')[1]}.repackaged.gz"
        output_coordinate_file_map[
            shard
        ] = f"{coordinates_output_dir}/{output_name}_coordinates.tsv"
    output_coordinate_files_exist = parallel_file_exists(
        list(output_coordinate_file_map.values())
    )

    shards_to_repackage = []
    for shard in output_coordinate_file_map:
        if not output_coordinate_files_exist[output_coordinate_file_map[shard]]:
            shards_to_repackage.append(shard)
    logger.info("Found %i shards to repackage", len(shards_to_repackage))

    # Process shards
    with run_batch(args, batch_name="Repackage VCF shards") as batch:
        for shard in tqdm(shards_to_repackage, unit="shards"):
            # Naming the output shard with this format:
            # ukb<field_id>_chr<#>_block<#>.repackaged.gz
            # Using 'repackaged' + gz extension since that is the extension DNANexus uses in their code
            output_name = f"ukb26041_{contig}_block{os.path.split(shard)[-1].split('.')[0].split('-')[1]}.repackaged.gz"
            output_tabix = f"{output_name}.tbi"
            output_coordinates = f"{output_name}_coordinates.tsv"
            j = init_job(
                batch,
                f"{shard} repackage",
                DOCKER_IMAGE,
                args.cpu,
                args.memory,
                args.disk_size,
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
                f"""bash ukbb_header_reformat.sh {local_vcf_path} {output_name} > {output_coordinates}"""
            )
            j.command(f"""gsutil -u maclab-ukbb -m cp {output_name} {output_dir}""")
            j.command(f"""gsutil -u maclab-ukbb -m cp {output_tabix} {output_dir}""")
            j.command(f"""gsutil -m cp {output_coordinates} {coordinates_output_dir}""")


if __name__ == "__main__":
    p = init_arg_parser(
        default_cpu=2, default_memory=8, default_billing_project="gnomad-production"
    )
    p.add_argument(
        "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    p.add_argument(
        "--contig", help="Chromosome to repackage",
    )
    p.add_argument(
        "--disk-size", help="Storage (in GB) for batch jobs", default=10, type=int
    )
    p.add_argument(
        "--output-dir",
        help="Where to write repackaged vcf output.",
        default="gs://broad-ukbb-requester-pays/broad.freeze_7/sharded_vcf/",
    )
    args = p.parse_args()
    main(args)
