import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    get_checkpoint_path,
    last_END_positions_ht_path,
    logging_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("rekey_last_end_ht")
logger.setLevel(logging.INFO)


def main(args):
    """
	This script is necessary because sparse MatrixTables now come keyed by locus and alleles by default.
	That means that `compute_last_ref_block_end` wrote our current last END HT keyed by locus and alleles.
    `compute_last_ref_block_end` (fixed with https://github.com/broadinstitute/gnomad_methods/pull/279) 
    now explicitly keys the last END HT by locus, so this script will not be necessary for any future 
    production projects.

	The last END HT is used to optimize densifies that only require densifying to a specific list of sites
	(as opposed to running a full `hl.experimental.densify` on all sites). 

	`densify_sites` annotates a HT with the desired sites with the last END HT using this code:
	`last_END_positions_ht[sites_ht.key].last_END_position`.

	Thus, we want the last END HT to be keyed by locus and not locus and alleles in case the sites HT does not
	contain alleles, or it contains a different set of alleles for one or more loci than the last END HT.
	"""
    hl.init(log="/rekey_last_end_ht.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze

    try:
        logger.info("Reading in last END HT and re-keying by locus...")
        ht = hl.read_table(last_END_positions_ht_path(freeze))
        ht = ht.key_by("locus")
        ht = ht.checkpoint(
            get_checkpoint_path(data_source, freeze, name="temp_last_END"),
            overwrite=True,
        )
        ht = ht.write(last_END_positions_ht_path(freeze), overwrite=True)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
