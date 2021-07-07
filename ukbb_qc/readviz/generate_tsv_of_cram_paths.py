import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    cram_map_path,
    readviz_haplotype_caller_path,
    readviz_ht_exploded_path,
    readviz_per_sample_tsv_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("get_tsv_for_haplotype_caller")
logger.setLevel(logging.INFO)


def main(args):
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    ht = hl.read_table(readviz_ht_exploded_path()).select()
    ht = ht.annotate(
        variants_tsv_bgz=hl.format("%s/%s.tsv.bgz", readviz_per_sample_tsv_path(), ht.S)
    )

    cram_ht = hl.import_table(cram_map_path(*tranche_data), no_header=True).key_by("f0")
    cram_ht = cram_ht.annotate(crai_path=hl.format("%s.crai", cram_ht.f1))
    cram_ht = cram_ht.rename({"f0": "S", "f1": "cram_path"})

    ht = ht.annotate(**cram_ht[ht.key])
    if ht.aggregate(hl.agg.count_where(hl.is_missing(ht.cram_path) > 0)):
        raise DataException(
            "Some samples are missing a cram path. Please double check and restart!"
        )
    if ht.aggregate(hl.agg.count_where(hl.is_missing(ht.variants_tsv_bgz) > 0)):
        raise DataException(
            "Some samples are missing their variants TSV files. Please double check and restart!"
        )

    ht = ht.rename({"S": "sample_id"})
    ht.export(
        f"{readviz_haplotype_caller_path()}/inputs/step4_output_cram_and_tsv_paths_table.tsv"
    )


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
