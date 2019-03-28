# Resources for UKBB QC

ukbb_array_dir = 'gs://broad-ukbb/data/array/'
ukbb_calling_intervals = 'gs://broad-ukbb/data/white_album_exome_calling_regions.v1.interval_list'
ukbb_calling_intervals_summary = 'gs://broad-ukbb/data/ukbb_exome_calling_intervals.summary.txt'
gc_content_bed = 'gs://broad-ukbb/data/gc5Base.bed'


def intervals_ht(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/coverage/intervals.ht'


def coverage_by_target(data_source: str, freeze: int = CURRENT_FREEZE, ht: bool = True) -> str:
    prefix = f'{sample_qc_prefix(data_source, freeze)}/coverage/coverage_by_target'
    if ht:
        return prefix + '.ht'
    else:
        return prefix + '.txt'

def coverage_by_region_type(data_source: str, freeze: int = CURRENT_FREEZE, gc_bin_type: str, ht: bool = True) -> str:
    prefix = f'{sample_qc_prefix(data_source, freeze)}/coverage/coverage_by_region_type.{gc_bin_type}'
    if ht:
        return prefix + '.ht'
    else:
        return prefix + '.txt'
