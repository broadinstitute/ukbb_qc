import argparse
import json
import logging

import hail as hl

from ukbb_qc.resources.basics import get_release_path, release_mt_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("get_vcf_shard_loci")
logger.setLevel(logging.INFO)


def main(args):
    """
    Get VCF shard position start and end positions.

    Use MatrixTable (MT) metadata to get first and last positions in each MT partition.
    Assumes MT was written with the same number of partitions as number of desired VCF shards.
    """

    hl.init(log="get_vcf_shard_loci.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    logger.info("Opening MatrixTable metadata file...")
    metadata_path = f"{release_mt_path(*tranche_data)}/rows/rows/metadata.json.gz"
    with hl.hadoop_open(metadata_path, "rb") as m:
        data = json.loads(m.read())
    partition_data = data["_jRangeBounds"]

    logger.info("Writing to output file...")
    shard_name_length = args.length
    output_text_file = f"{get_release_path(*tranche_data)}/shard_positions.tsv"
    with hl.hadoop_open(output_text_file, "w") as o:
        o.write("shard_name\tstart\tend\n")
        for index, part_data in enumerate(partition_data):
            shard_int = format(index, f"0{shard_name_length}")
            shard_name = f"part-{shard_int}.bgz"
            start_locus = f"{part_data['start']['locus']['contig']}:{part_data['start']['locus']['position']}"
            end_locus = f"{part_data['end']['locus']['contig']}:{part_data['end']['locus']['position']}"
            o.write(f"{shard_name}\t{start_locus}\t{end_locus}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script searches for regional missense constraint in gnomAD"
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    args = parser.parse_args()
    main(args)
