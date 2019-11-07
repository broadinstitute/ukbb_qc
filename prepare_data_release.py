from gnomad_hail import *
from gnomad_hail.resources.variant_qc import *
from gnomad_qc.variant_qc.prepare_data_release import * 
from ukbb_qc.resources import *
import copy
import itertools


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("data_release")
logger.setLevel(logging.INFO)


GROUPS = ['adj', 'raw']
SEXES = ['male', 'female', 'ambiguous_sex', 'sex_aneuploidy']
POPS = ['afr', 'amr', 'asj', 'eas', 'fin', 'nfe', 'oth', 'sas']
FAF_POPS = ['afr', 'amr', 'eas', 'nfe', 'sas']
# hybrid pops (freeze 4): 0 (afr -- Caribbean/African), 1 (sas/eas -- Chinese), 2 (afr -- White and Black African/Caribbean), 3 (sas/eas -- Indian/Pakistani)
# 4: (asj/nfe/amr -- British/other white), 5: (nfe/amr --other/other Asian), 6 (nfe -- British), 7 (nfe/fin -- British), 8 (nfe/fin -- other white)
# 9: (nfe -- British/Irish)
NFE_SUBPOPS = ['4', '6', '7', '8', '9']
AFR_SUBPOPS = ['0', '2']
SAS_SUBPOPS = ['3']
EAS_SUBPOPS = ['1']
SORT_ORDER = ['popmax', 'group', 'pop', 'subpop', 'sex']


def sample_sum_check(ht: hl.Table, prefix: str, label_groups: Dict[str, List[str]], verbose: bool, subpop=None):
    '''
    Compute afresh the sum of annotations for a specified group of annotations, and compare to the annotated version;
    display results from checking the sum of the specified annotations in the terminal
    :param Table ht: Hail Table containing annotations to be summed
    :param str prefix: Subset of gnomAD
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks
    :param str subpop: Subpop abbreviation, supplied only if subpopulations are included in the annotation groups being checked
    :rtype: None
    '''
    combo_AC = [ht.info[f'{prefix}_AC_{x}'] for x in make_label_combos(label_groups)]
    combo_AN = [ht.info[f'{prefix}_AN_{x}'] for x in make_label_combos(label_groups)]
    combo_nhomalt = [ht.info[f'{prefix}_nhomalt_{x}'] for x in make_label_combos(label_groups)]
    group = label_groups.pop('group')[0]
    alt_groups = "_".join(sorted(label_groups.keys(), key=lambda x: SORT_ORDER.index(x)))

    annot_dict = {f'sum_AC_{group}_{alt_groups}': hl.sum(combo_AC),
                  f'sum_AN_{group}_{alt_groups}': hl.sum(combo_AN),
                  f'sum_nhomalt_{group}_{alt_groups}': hl.sum(combo_nhomalt)}

    ht = ht.annotate(**annot_dict)

    for subfield in ['AC', 'AN', 'nhomalt']:
        if not subpop:
            generic_field_check(ht, (ht.info[f'{prefix}_{subfield}_{group}'] != ht[f'sum_{subfield}_{group}_{alt_groups}']),
                                f'{prefix}_{subfield}_{group} = sum({subfield}_{group}_{alt_groups})',
                                [f'info.{prefix}_{subfield}_{group}', f'sum_{subfield}_{group}_{alt_groups}'], verbose)
        else:
            generic_field_check(ht, (ht.info[f'{prefix}_{subfield}_{group}_{subpop}'] != ht[f'sum_{subfield}_{group}_{alt_groups}']),
                                f'{prefix}_{subfield}_{group}_{subpop} = sum({subfield}_{group}_{alt_groups})',
                                [f'info.{prefix}_{subfield}_{group}_{subpop}', f'sum_{subfield}_{group}_{alt_groups}'], verbose)


def generic_field_check(ht: hl.Table, cond_expr, check_description, display_fields, verbose):
    '''
    Check a generic logical condition involving annotations in a Hail Table and print the results to terminal
    :param Table ht: Table containing annotations to be checked
    :param Expression cond_expr: logical expression referring to annotations in ht to be checked
    :param str check_description: String describing the condition being checked; is displayed in terminal summary message
    :param list of str display_fields: List of names of ht annotations to be displayed in case of failure (for troubleshooting purposes);
        these fields are also displayed if verbose is True
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks
    :rtype: None
    '''
    ht_orig = ht
    ht = ht.filter(cond_expr)
    n_fail = ht.count()
    if n_fail > 0:
        logger.info(f'Found {n_fail} sites that fail {check_description} check:')
        ht = ht.flatten()
        ht.select('locus', 'alleles', *display_fields).show()
    else:
        logger.info(f'PASSED {check_description} check')
        if verbose:
            ht_orig = ht_orig.flatten()
            ht_orig.select(*display_fields).show()


def sanity_check_ht(ht: hl.Table, missingness_threshold=0.5, verbose=False):
    '''
    Perform a battery of sanity checks on a specified group of subsets in a Hail Table containing variant annotations;
    includes summaries of % filter status for different partitions of variants; histogram outlier bin checks; checks on
    AC, AN, and AF annotations; checks that subgroup annotation values add up to the supergroup annotation values;
    checks on sex-chromosome annotations; and summaries of % missingness in variant annotations

    :param Table ht: Table containing variant annotations to check
    :param bool verbose: If True, display top values of relevant annotations being checked, regardless of whether check
        conditions are violated; if False, display only top values of relevant annotations if check conditions are violated
    :return: Terminal display of results from the battery of sanity checks
    :rtype: None
    '''
    n_sites = ht.count()
    contigs = ht.aggregate(hl.agg.collect_as_set(ht.locus.contig))
    logger.info(f'Found {n_sites} sites in {data_source} for contigs {contigs}')
    info_metrics = list(ht.row.info)
    non_info_metrics = list(ht.row)
    non_info_metrics.remove('info')

    logger.info('VARIANT FILTER SUMMARIES:')
    ht = ht.annotate(is_filtered=ht.filters.length() > 0,
                     in_problematic_region=hl.any(lambda x: x, [ht.info.lcr, ht.info.segdup, ht.info.decoy]))

    ht_filter_check1 = ht.group_by(ht.is_filtered).aggregate(**make_filters_sanity_check_expr(ht)).order_by(hl.desc('n'))
    ht_filter_check1.show()

    ht_filter_check2 = ht.group_by(ht.info.allele_type).aggregate(**make_filters_sanity_check_expr(ht)).order_by(hl.desc('n'))
    ht_filter_check2.show()

    ht_filter_check3 = (ht.group_by(ht.info.allele_type, ht.in_problematic_region)
                        .aggregate(**make_filters_sanity_check_expr(ht)).order_by(hl.desc('n')))
    ht_filter_check3.show(50, 140)

    ht_filter_check4 = (ht.group_by(ht.info.allele_type, ht.in_problematic_region, ht.info.n_alt_alleles)
                        .aggregate(**make_filters_sanity_check_expr(ht)).order_by(hl.desc('n')))
    ht_filter_check4.show(50, 140)

    logger.info('HISTOGRAM CHECKS:')
    for hist in ['gq_hist_alt', 'gq_hist_all', 'dp_hist_alt', 'dp_hist_all', 'ab_hist_alt']:
        # Check subfield == 0
        generic_field_check(ht, (ht.info[f'{hist}_n_smaller'] != 0), f'{hist}_n_smaller == 0',
                            [f'info.{hist}_n_smaller'], verbose)
        if hist not in ['dp_hist_alt', 'dp_hist_all']:  # NOTE: DP hists can have nonzero values in n_larger bin
            generic_field_check(ht, (ht.info[f'{hist}_n_larger'] != 0), f'{hist}_n_larger == 0',
                                [f'info.{hist}_n_larger'], verbose)

    logger.info('RAW AND ADJ CHECKS:')
    for subfield in ['AC', 'AN', 'AF']:
        # Check AC, AN, AF > 0
        generic_field_check(ht, (ht.info[f'{subfield}_raw'] <= 0), f'{subfield}_raw > 0',
                            [f'info.{subfield}_raw'], verbose)
        generic_field_check(ht, (ht.info[f'{subfield}_adj'] < 0), f'{subfield}_adj >= 0',
                            [f'info.{subfield}_adj', 'filters'], verbose)

    for subfield in ['AC', 'AN', 'nhomalt']:
        # Check AC_raw >= AC adj
        generic_field_check(ht, (ht.info[f'{subfield}_raw'] < ht.info[f'{subfield}_adj']),
                            f'{subfield}_raw >= {subfield}_adj',
                            [f'info.{subfield}_raw', f'info.{subfield}_adj'], verbose)

    logger.info('SAMPLE SUM CHECKS:')
    # NOTE: group entry in dict is required here to be a list of length 1
    sample_sum_check(ht, dict(group=['adj'], pop=pop_adjusted), verbose)
    sample_sum_check(ht, dict(group=['adj'], sex=SEXES), verbose)
    sample_sum_check(ht, dict(group=['adj'], pop=pop_adjusted, sex=SEXES), verbose)

    logger.info('SEX CHROMOSOME ANNOTATION CHECKS:')
    female_metrics = [x for x in info_metrics if '_female' in x]
    male_metrics = [x for x in info_metrics if '_male' in x]

    if 'Y' in contigs:
        logger.info('Check values of female metrics for Y variants are NA:')
        ht_y = hl.filter_intervals(ht, [hl.parse_locus_interval('Y')])
        metrics_values = {}
        for metric in female_metrics:
            metrics_values[metric] = hl.agg.collect_as_set(ht_y.info[metric])
        output = ht_y.aggregate(hl.struct(**metrics_values))
        for metric,values in dict(output).items():
            if values == {None}:
                logger.info(f"PASSED {metric} = {None} check for Y variants")
            else:
                logger.info(f"FAILED Y check: Found {values} in {metric}")

    logger.info('Check values of male nhomalt metrics for X nonpar variants are 0:')
    ht_x = hl.filter_intervals(ht, [hl.parse_locus_interval('X')])
    ht_xnonpar = ht_x.filter(ht_x.locus.in_x_nonpar())
    n = ht_xnonpar.count()
    logger.info(f"Found {n} X nonpar sites")  # Lots of values found in male X nonpar sites

    male_metrics = [x for x in male_metrics if 'nhomalt' in x]
    metrics_values = {}
    for metric in male_metrics:
        metrics_values[metric] = hl.agg.collect_as_set(ht_xnonpar.info[metric])
    output = ht_xnonpar.aggregate(hl.struct(**metrics_values))
    for metric, values in dict(output).items():
        if values == {0}:
            logger.info(f"PASSED {metric} = 0 check for X nonpar variants")
        else:
            logger.info(f"FAILED X nonpar check: Found {values} in {metric}")

    logger.info('Check (nhomalt == nhomalt_female) for X nonpar variants:')
    female_metrics = [x for x in female_metrics if 'nhomalt' in x]
    for metric in female_metrics:
        standard_field = metric.replace('_female', '')
        generic_field_check(ht_xnonpar, (ht_xnonpar.info[f'{metric}'] != ht_xnonpar.info[f'{standard_field}']),
                            f'{metric} == {standard_field}', [f'info.{metric}', f'info.{standard_field}'], verbose)

    logger.info('MISSINGNESS CHECKS:')
    metrics_frac_missing = {}
    for x in info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht.info[x])) / n_sites
    for x in non_info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht[x])) / n_sites
    output = ht.aggregate(hl.struct(**metrics_frac_missing))

    n_fail = 0
    for metric, value in dict(output).items():
        if value > missingness_threshold:
            logger.info("FAILED missing check for {}: {}% missing".format(metric, 100 * value))
            n_fail += 1
        else:
            logger.info("Missingness check for {}: {}% missing".format(metric, 100 * value))
    logger.info("{} missing metrics checks failed".format(n_fail))


def make_freq_meta_index_dict(freq_meta):
    '''
    Make dictionary of the entries in the frequency array annotation, where keys are the grouping combinations and the values
    are the 0-based integer indices
    :param list of str freq_meta: Ordered list containing string entries describing all the grouping combinations contained in the
        frequency array annotation
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    '''
    index_dict = index_globals(freq_meta, dict(group=GROUPS))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=POPS)))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, sex=SEXES)))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=POPS, sex=SEXES)))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=['nfe'], subpop=NFE_SUBPOPS)))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=['eas'], subpop=EAS_SUBPOPS)))
    return index_dict


def make_info_dict(label_groups=None, bin_edges=None, faf=False, popmax=False, age_hist_data=None):
    '''
    Generate dictionary of Number and Description attributes to be used in the VCF header
    
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :param dict bin_edges: Dictionary keyed by annotation type, with values that reflect the bin edges corresponding to the annotation
    :param bool faf: If True, use alternate logic to auto-populate dictionary values associated with filter allele frequency annotations
    :param bool popmax: If True, use alternate logic to auto-populate dictionary values associated with popmax annotations
    :param str age_hist_data: Pipe-delimited string of age histograms, from get_age_distributions (somewhat required if popmax == True)
    :return: Dictionary keyed by VCF INFO annotations, where values are Dictionaries of Number and Description attributes
    :rtype: Dict of str: (Dict of str: str)
    '''
    info_dict = dict()

    if popmax:
        popmax_text = ""
        popmax_dict = {
            'popmax': {"Number": "A",
                                 "Description": "Population with maximum AF{}".format(popmax_text)},
            'AC_popmax': {"Number": "A",
                                    "Description": "Allele count in the population with the maximum AF{}".format(popmax_text)},
            'AN_popmax': {"Number": "A",
                                    "Description": "Total number of alleles in the population with the maximum AF{}".format(popmax_text)},
            'AF_popmax': {"Number": "A",
                                    "Description": "Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry){}".format(popmax_text)},
            'nhomalt_popmax': {"Number": "A",
                                         "Description": "Count of homozygous individuals in the population with the maximum allele frequency{}".format(popmax_text)}
        }
        info_dict.update(popmax_dict)
        age_hist_dict = {
            "age_hist_het_bin_freq": {"Number": "A",
                                                "Description": f"Histogram of ages of heterozygous individuals; bin edges are: {bin_edges[f'{prefix}_het']}; total number of individuals of any genotype bin: {age_hist_data}"},
            "age_hist_het_n_smaller": {"Number": "A",
                                                 "Description": "Count of age values falling below lowest histogram bin edge for heterozygous individuals"},
            "age_hist_het_n_larger": {"Number": "A",
                                                "Description": "Count of age values falling above highest histogram bin edge for heterozygous individuals"},
            "age_hist_hom_bin_freq": {"Number": "A",
                                                "Description": f"Histogram of ages of homozygous alternate individuals; bin edges are: {bin_edges[f'{prefix}_hom']}; total number of individuals of any genotype bin: {age_hist_data}"},
            "age_hist_hom_n_smaller": {"Number": "A",
                                                 "Description": "Count of age values falling below lowest histogram bin edge for homozygous alternate individuals"},
            "age_hist_hom_n_larger": {"Number": "A",
                                                "Description": "Count of age values falling above highest histogram bin edge for homozygous alternate individuals"}
        }
        info_dict.update(age_hist_dict)
    return info_dict


def make_hist_bin_edges_expr(ht):
    '''
    Create dictionary containing variant histogram annotations and their associated bin edges, formatted into a string
    separated by pipe delimiters

    :param Table ht: Table containing histogram variant annotations
    :return: Dictionary keyed by histogram annotation name, with corresponding reformatted bin edges for values
    :rtype: Dict of str: str
    '''
    edges_dict = {'het': '|'.join(map(lambda x: f'{x:.1f}', ht.take(1)[0].age_hist_het[0].bin_edges)),
                  'hom': '|'.join(map(lambda x: f'{x:.1f}', ht.take(1)[0].age_hist_hom[0].bin_edges))}
    for hist in HISTS:
        edges_dict[hist] = '|'.join(map(lambda x: f'{x:.2f}', ht.take(1)[0][hist].bin_edges)) if 'ab' in hist else \
            '|'.join(map(lambda x: str(int(x)), ht.take(1)[0][hist].bin_edges))
    return edges_dict


def make_index_dict(ht):
    '''
    Create a look-up Dictionary for entries contained in the frequency annotation array
    :param Table ht: Table containing freq_meta global annotation to be indexed
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    '''
    freq_meta = hl.eval(ht.globals.freq_meta)
    index_dict = make_freq_meta_index_dict(freq_meta)
    return index_dict


def get_age_distributions(data_source, freeze):
    """
    Get background distribution of age (field 21022 is age at recruitment)

    :param str data_source: One of regeneron or broad
    :param int freeze: One of data freezes
    :return: pipe-delimited string with ages in pre-determined bins (<30, 30-35, ..., 75-80, 80+)
    :rtype: str
    """
    ukbb_phenotypes = hl.import_table(ukbb_phenotype_path, impute=True)
    ukbb_phenotypes = ukbb_phenotypes.key_by(s_old=hl.str(ukbb_phenotypes['f.eid']))
    ukbb_age = ukbb_phenotypes.select('f.21022.0.0')
    sample_map_ht = hl.read_table(array_sample_map_ht(data_source, freeze))
    sample_map_ht = sample_map_ht.key_by('ukbb_app_26041_id')
    ukbb_age = ukbb_age.key_by(s=sample_map_ht[ukbb_age.key].s)

    age_hist_data = ukbb_age.aggregate(hl.agg.hist(ukbb_age['f.21022.0.0'], 30, 80, 10))
    age_hist_data.bin_freq.insert(0, age_hist_data.n_smaller)
    age_hist_data.bin_freq.append(age_hist_data.n_larger)
    return age_hist_data.bin_freq


def main(args):
    hl.init(log='/release.log', default_reference='GRCh38')

    data_source = args.data_source
    freeze = args.freeze
    age_hist_data = get_age_distributions(data_source, freeze)

    if args.prepare_internal_ht:
        
        freq_ht = hl.read_table(var_annotations_ht_path(data_source, 'join_freq'))
        index_dict = make_index_dict(freq_ht)
        rf_ht = hl.read_table(var_annotations_ht_path(data_source, 'rf')).drop('info_ac', 'ac', 'ac_raw')
        vep_ht = hl.read_table(var_annotations_ht_path(data_source, 'vep'))
        dbsnp_ht = hl.read_table(dbsnp_ht_path)
        hist_ht = hl.read_table(var_annotations_ht_path(data_source, 'qual_hists'))
        hist_ht = hist_ht.select('gq_hist_alt', 'gq_hist_all', 'dp_hist_alt', 'dp_hist_all', 'ab_hist_alt')
        allele_ht = hl.read_table(var_annotations_ht_path(data_source, 'allele_data'))

        logger.info('Adding annotations...')
        ht = prepare_table_annotations(freq_ht, rf_ht, vep_ht, dbsnp_ht, hist_ht, index_dict, allele_ht)
        ht = ht.annotate_globals(freq_index_dict=make_index_dict(ht), faf_index_dict=make_faf_index_dict(ht),
                                 age_distribution=age_hist_data)
        ht.write(release_ht_path(data_source, freeze), args.overwrite)


    if args.prepare_release_vcf:
        ht = hl.read_table(release_ht_path(data_source, freeze))
        bin_edges = make_hist_bin_edges_expr(ht)

        # Make INFO dictionary for VCF
        INFO_DICT.update(make_info_dict(bin_edges=bin_edges, popmax=True,
                                        age_hist_data='|'.join(str(x) for x in age_hist_data)))
        INFO_DICT.update(make_info_dict(dict(group=GROUPS)))
        INFO_DICT.update(make_info_dict(dict(group=GROUPS, sex=SEXES)))
        INFO_DICT.update(make_info_dict(dict(group=GROUPS, pop=POPS)))
        INFO_DICT.update(make_info_dict(dict(group=GROUPS, pop=POPS, sex=SEXES)))
        INFO_DICT.update(make_info_dict(dict(group=GROUPS, pop=['nfe'], subpop=NFE_SUBPOPS)))
        INFO_DICT.update(make_info_dict(dict(group=GROUPS, pop=['eas'], subpop=EAS_SUBPOPS)))
        INFO_DICT.update(make_info_dict(dict(group=['adj']), faf=True))
        INFO_DICT.update(make_info_dict(dict(group=['adj'], pop=FAF_POPS), faf=True))
        INFO_DICT.update(make_hist_dict(bin_edges))

        # Adjust keys to remove gnomad, adj tags before exporting to VCF
        #new_info_dict = {i.replace('gnomad_', '').replace('_adj', ''): j for i,j in INFO_DICT.items()}
        new_info_dct = {i.replace('_adj', ''): j for i,j in INFO_DICT.items()}

        # Construct INFO field
        ht = ht.annotate(info=hl.struct(**make_info_expr(ht)))
        ht = ht.annotate(info=ht.info.annotate(**unfurl_nested_annotations(ht)))
        ht = set_female_y_metrics_to_na(ht)

        # Select relevant fields for VCF export
        ht = ht.select('info', 'filters', 'rsid', 'qual', 'vep')
        ht.write(release_ht_path(data_source, nested=False, temp=True), args.overwrite)

        # Move 'info' annotations to top level for browser release
        ht = hl.read_table(release_ht_path(data_source, nested=False, temp=True))
        ht = ht.transmute(**ht.info)
        ht = ht.select_globals('rf')
        ht.write(release_ht_path(data_source, nested=False), args.overwrite)

        # Remove gnomad_ prefix for VCF export
        ht = hl.read_table(release_ht_path(data_source, nested=False, temp=True))
        rg = hl.get_reference('GRCh38')
        contigs = rg.contigs[:24] # autosomes + X/Y
        #contigs = hl.eval(ht.aggregate(hl.agg.collect_as_set(ht.locus.contig)))
        ht = ht.drop('vep')
        row_annots = list(ht.row.info)
        #new_row_annots = [x.replace('gnomad_', '').replace('_adj', '') for x in row_annots]
        new_row_annots = [x.replace('_adj', '') for x in row_annots]
        info_annot_mapping = dict(zip(new_row_annots, [ht.info[f'{x}'] for x in row_annots]))
        ht = ht.transmute(info=hl.struct(**info_annot_mapping))

        # Rearrange INFO field in desired ordering
        drop_hists = [x + '_n_smaller' for x in HISTS] + [x + '_bin_edges' for x in HISTS] + [x + '_n_larger' for x in HISTS if 'dp_' not in x] + ['age_hist_hom_bin_edges', 'age_hist_het_bin_edges']
        ht = ht.annotate(info=ht.info.select('AC', 'AN', 'AF', 'rf_tp_probability',
                                             *ht.info.drop('AC', 'AN', 'AF', 'rf_tp_probability', *drop_hists)))

        # Add VEP annotations
        vep_csq_ht = hl.read_table(annotations_ht_path(data_source, 'vep_csq'))
        new_info_dict.update({'vep': {'Description': hl.eval(vep_csq_ht.globals.vep_csq_header)}})
        header_dict = {'info': new_info_dict,
                       'filter': make_filter_dict(ht)}
        ht = ht.annotate(info=ht.info.annotate(vep=vep_csq_ht[ht.key].vep))

        # Export VCFs, full and by chromosome
        #gnomad_ref = hl.ReferenceGenome.read('gs://gnomad-public/resources/gnomad_grch37.json')
        ht = ht.key_by(locus=hl.locus(ht.locus.contig, ht.locus.position),
                       alleles=ht.alleles)

        for contig in contigs:
            contig_ht = hl.filter_intervals(ht, [hl.parse_locus_interval(contig)])
            mt = hl.MatrixTable.from_rows_table(contig_ht).key_cols_by(s='foo')
            hl.export_vcf(mt, release_vcf_path(data_source, contig=contig), metadata=header_dict)

        mt = hl.MatrixTable.from_rows_table(ht).key_cols_by(s='foo')
        hl.export_vcf(mt, release_vcf_path(data_source), metadata=header_dict)

    
    if args.sanity_check_sites:
        ht = hl.read_table(release_ht_path(data_source, freeze, nested=False, temp=True))
        sanity_check_ht(ht, missingness_threshold=0.5, verbose=args.verbose)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    parser.add_argument('--prepare_internal_ht', help='Prepare internal HailTable', action='store_true')
    parser.add_argument('--prepare_release_vcf', help='Prepare release VCF', action='store_true')
    parser.add_argument('--sanity_check_sites', help='Run sanity checks function', action='store_true')
    parser.add_argument('--verbose', help='Run sanity checks function with verbose output', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
