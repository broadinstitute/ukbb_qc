import asyncio
import argparse
import logging
import os
from typing import Dict, List
from concurrent.futures import ThreadPoolExecutor

import hailtop.batch as hb
from hailtop.aiotools import LocalAsyncFS, RouterAsyncFS, AsyncFS
from hailtop.aiogoogle import GoogleStorageAsyncFS
from hailtop.utils import bounded_gather, tqdm
import hail as hl

from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    readviz_ht_exploded_path,
    readviz_per_sample_tsv_path,
)
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("run_per_sample_tsv_export")
logger.setLevel(logging.INFO)


async def parallel_file_exists(fnames: List[str]) -> Dict[str, bool]:
    """
    Check if file exists.

    Code is from Dan King to speed up check for files 
    (`file_exists` in gnomAD repo was too slow for large number of Batch jobs).
    """

    async def low_level_async_file_exists(fs: AsyncFS, url: str):
        try:
            await fs.statfile(url)
        except FileNotFoundError:
            return False
        else:
            return True

    async def async_file_exists(fs: AsyncFS, fname: str):
        fext = os.path.splitext(fname)[1]
        if fext in [".ht", ".mt"]:
            fname += "/_SUCCESS"
        return await low_level_async_file_exists(fs, fname)

    with tqdm(
        total=len(fnames), desc="check files for existence", disable=False
    ) as pbar:
        with ThreadPoolExecutor() as thread_pool:
            async with RouterAsyncFS(
                "file", [LocalAsyncFS(thread_pool), GoogleStorageAsyncFS()]
            ) as fs:

                def create_unapplied_function(fname):
                    async def unapplied_function():
                        x = await async_file_exists(fs, fname)
                        pbar.update(1)
                        return x

                    return unapplied_function

                file_existence_checks = [
                    create_unapplied_function(fname) for fname in fnames
                ]
                file_existence = await bounded_gather(
                    *file_existence_checks, parallelism=750
                )
    return {fname: exists for fname, exists in zip(fnames, file_existence)}


def get_sample_ids(ids_file: str, header: bool = False) -> List[str]:
    """
    Open file with all sample IDs and stores IDs in a list.

    :param str ids_file: Path to file containing sample IDs with variants to extract for readviz.
    :param bool header: Whether IDs file has a header line. Default is False.
    :return: List of sample IDs
    :rtype: List[str]
    """
    sample_ids = []
    with open(ids_file) as i:
        if header:
            header = i.readline()
        for line in i:
            sample_ids.append(line.strip())
    return sample_ids


def export_tsv(ht_path: str, sample_id: str, tsv_path: str, success_path: str) -> None:
    """
    Read in hail Table, filter to specified sample ID, and export TSV.

    Also write an empty success file to ensure no Batch jobs finish with only partially exported TSV.

    :param str ht_path: Path to input hail Table.
    :param str sample_id: Sample for which to export a TSV.
    :param str tsv_path: Path to output TSV.
    :param str success_path: Path to output success file.
    :return: None
    """
    ht = hl.read_table(ht_path)
    ht = ht.filter(ht.S == sample_id)
    ht = ht.naive_coalesce(1)
    ht.export(tsv_path)

    with hl.hadoop_open(success_path, "w") as o:
        o.write("")


def main(args):

    backend = hb.ServiceBackend("gnomad-production", "gnomad-kc")
    b = hb.Batch(
        backend=backend,
        default_cpu=1,
        default_python_image="gcr.io/broad-mpg-gnomad/tgg-methods-vm:20210623",
        project="broad-mpg-gnomad",
    )

    logger.info("Extracting sample IDs...")
    sample_ids = get_sample_ids(args.ids_file, args.header)

    logger.info("Preparing to start batch job...")
    files = [
        f"{readviz_per_sample_tsv_path()}/{sample}_success.txt" for sample in sample_ids
    ]
    file_exists = asyncio.get_event_loop().run_until_complete(
        parallel_file_exists(files)
    )

    for sample in sample_ids:
        logger.info("Working on %s", sample)
        if file_exists[f"{readviz_per_sample_tsv_path()}/{sample}_success.txt"]:
            logger.info(
                "Output success txt file already exists, skipping %s...", sample
            )
            continue
        j = b.new_python_job(name=sample)
        j.call(
            export_tsv,
            readviz_ht_exploded_path(),
            sample,
            f"{readviz_per_sample_tsv_path()}/{sample}.tsv",
            f"{readviz_per_sample_tsv_path()}/{sample}_success.txt",
        )

    b.run(wait=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids-file", help="File with sample IDs")
    parser.add_argument(
        "--header",
        help="Whether file with sample IDs has a header line",
        action="store_true",
    )
    parser.add_argument("--slack-channel", help="Send message to Slack channel/user")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
