import asyncio
import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import parallel_file_exists
from gnomad.utils.slack import slack_notifications


from ukbb_qc.utils.utils import get_sample_ids
from ukbb_qc.resources.basics import (
    cram_map_path,
    readviz_haplotype_caller_path,
    readviz_ids_tsv_path,
    readviz_per_sample_tsv_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("get_tsv_for_haplotype_caller")
logger.setLevel(logging.INFO)


def main(args):
    freeze = args.freeze

    logger.info("Getting sample IDs...")
    sample_ids = get_sample_ids(readviz_ids_tsv_path(freeze))

    logger.info("Starting TSV file existence checks...")
    tsvs = [
        f"{readviz_per_sample_tsv_path()}/{sample}.tsv.bgz" for sample in sample_ids
    ]
    tsv_file_exists = asyncio.get_event_loop().run_until_complete(
        parallel_file_exists(tsvs)
    )

    logger.info("Starting cram existence checks...")
    cram_map = {}
    with hl.hadoop_open(cram_map_path(freeze)) as c:
        for line in c:
            sample, cram = line.strip().split("\t")
            cram_map[sample] = cram
    crams_exist = asyncio.get_event_loop().run_until_complete(
        parallel_file_exists(list(cram_map.values()))
    )

    logger.info("Starting to write to output TSV...")
    with hl.hadoop_open(readviz_ids_tsv_path(freeze)) as s, hl.hadoop_open(
        f"{readviz_haplotype_caller_path()}/inputs/step4_output_cram_and_tsv_paths_table.tsv",
        "w",
    ) as o:
        o.write("sample_id\tcram\tcrai\tvariants_tsv_bgz\n")

        for line in s:
            sample = line.strip()
            if sample not in cram_map:
                raise DataException(
                    f"{sample} is missing a cram path. Please double check and restart!"
                )
            else:
                cram = cram_map[sample]
                tsv = f"{readviz_per_sample_tsv_path()}/{sample}.tsv.bgz"
                if not crams_exist[cram]:
                    raise DataException(
                        f"{sample}'s cram does not exist. Please double check and restart!"
                    )
            if not tsv_file_exists[tsv]:
                raise DataException(
                    f"{sample} is missing their variants TSV file. Please double check and restart!"
                )
            o.write(f"{sample}\t{cram}\t{cram}.crai\t{tsv}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument("--slack_channel", help="Send message to Slack channel/user")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
