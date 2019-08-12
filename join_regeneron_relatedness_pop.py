from gnomad_hail import *
from gnomad_hail.utils.sample_qc import *
from ukbb_qc.resources import *
import hail as hl
import argparse

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("join_regeneron_relatedness_pop")
logger.setLevel(logging.INFO)


def main(args):

    freeze = args.freeze
    data_source = args.data_source

    if args.join_relatedness:
        # Define Relationship thresholds
        first_degree_threshold = [float(i) for i in args.first_degree_kin_threshold.split(",")]
        second_degree_threshold = args.second_degree_kin_threshold
        ibd2_parent_offspring_threshold = args.ibd2_parent_offspring_threshold

        # Load Regeneron relationship info
        regeneron_rel_SD_ht = hl.import_table(get_regeneron_relatedness_path(freeze, relationship='2nd-degree'), impute=True, delimiter=" ")
        regeneron_rel_FS_ht = hl.import_table(get_regeneron_relatedness_path(freeze, relationship='full-sibling'), impute=True, delimiter=" ")
        regeneron_rel_PC_ht = hl.import_table(get_regeneron_relatedness_path(freeze, relationship='parent-child'), impute=True, delimiter=" ")
        regeneron_rel_ht = regeneron_rel_SD_ht.union(regeneron_rel_FS_ht).union(regeneron_rel_PC_ht)

        # Add a new column to each of the tables that sorts and joins the ids to be used to index pairs
        relatedness_ht = hl.read_table(relatedness_ht_path(data_source, freeze))
        relatedness_ht = relatedness_ht.annotate(samples=hl.delimit(hl.sorted([relatedness_ht.i.s, relatedness_ht.j.s])))
        regeneron_rel_ht = regeneron_rel_ht.annotate(samples=hl.delimit(hl.sorted([regeneron_rel_ht.IID1, regeneron_rel_ht.IID2])))

        # Add relationship status similar to Regeneron's based on thresholds above and join with Regeneron table
        relatedness_ht = relatedness_ht.annotate(Broad_classification=hl.case()
                                                 .when((relatedness_ht.kin > second_degree_threshold) &
                                                       (relatedness_ht.kin < first_degree_threshold[0]), 'Second-degree')
                                                 .when((relatedness_ht.kin > first_degree_threshold[0]) &
                                                       (relatedness_ht.kin < first_degree_threshold[1]) &
                                                       (relatedness_ht.ibd2 >= ibd2_parent_offspring_threshold), 'Full-sibling')
                                                 .when((relatedness_ht.kin > first_degree_threshold[0]) &
                                                       (relatedness_ht.kin < first_degree_threshold[1]) &
                                                       (relatedness_ht.ibd2 < ibd2_parent_offspring_threshold), 'Parent-child')
                                                 .default("None"),
                                                 Broad_IBD0=relatedness_ht.ibd0,
                                                 Broad_IBD1=relatedness_ht.ibd1,
                                                 Broad_IBD2=relatedness_ht.ibd2,
                                                 Broad_kin=relatedness_ht.kin)

        regeneron_rel_ht = regeneron_rel_ht.annotate(Regeneron_classification=(hl.switch(regeneron_rel_ht.MOST_LIKELY_REL)
                                                                               .when("CGH", 'Second-degree')
                                                                               .when("HAG", 'Second-degree')
                                                                               .when("FS", 'Full-sibling')
                                                                               .when("PC", 'Parent-child')
                                                                               .default("None")),
                                                     Regeneron_IBD0=regeneron_rel_ht.IBD0,
                                                     Regeneron_IBD1=regeneron_rel_ht.IBD1,
                                                     Regeneron_IBD2=regeneron_rel_ht.IBD2)

        relatedness_joined_ht = relatedness_ht.key_by('samples').join(regeneron_rel_ht.key_by('samples'), 'outer')
        relatedness_joined_ht = relatedness_joined_ht.annotate(Regeneron_classification=hl.or_else(relatedness_joined_ht.Regeneron_classification, "None"),
                                                               Broad_classification=hl.or_else(relatedness_joined_ht.Broad_classification, "None"))

        relatedness_joined_ht.write(get_regeneron_broad_relatedness_path(data_source, freeze), overwrite=args.overwrite)


    if args.filter_ukbb_ancestry:
        ukbb_phenotypes = hl.import_table(ukbb_phenotype_path , impute=True)
        ukbb_phenotypes = ukbb_phenotypes.key_by(s_old=hl.str(ukbb_phenotypes['f.eid']))
        ukbb_ancestry = ukbb_phenotypes.select('f.21000.0.0', 'f.21000.1.0', 'f.21000.2.0')
        sample_map = hl.import_table(array_sample_map, delimiter=',')
        sample_map = sample_map.key_by(s=hl.str(sample_map.eid_26041))
        ukbb_ancestry = ukbb_ancestry.key_by(s=sample_map[ukbb_ancestry.key].eid_sample)
        ukbb_ancestry = ukbb_ancestry.annotate(Self_reported_ancestry=hl.switch(ukbb_ancestry['f.21000.0.0'])
                                       .when(1,"White").when(1001,"British").when(1002,"Irish").when(1003,"Other white")
                                       .when(2,"Mixed").when(2001,"White and Black Caribbean").when(2002,"White and Black African").when(2003,"White and Asian")
                                       .when(2004,"Other Mixed").when(3,"Asian or Asian British").when(3001,"Indian").when(3002,"Pakistani")
                                       .when(3003,"Bangladeshi").when(3004,"Other Asian").when(4,"Black or Black British").when(4001,"Caribbean")
                                       .when(4002,"African").when(4003,"Other Black").when(5,"Chinese").when(6,"Other").when(-1,"Do not Know")
                                       .when(-3,"Prefer not to answer").default("None"))
        ukbb_ancestry.write(get_ukbb_self_reported_ancestry_path(freeze), overwrite=args.overwrite)


    if args.join_ancestry:
        ukbb_ancestry = hl.read_table(get_ukbb_self_reported_ancestry_path(freeze))
        regeneron_pops_ht = hl.import_table(get_regeneron_ancestry_path(freeze), impute=True).key_by('IID')
        regeneron_pops_ht = regeneron_pops_ht.annotate(Regeneron_pop=regeneron_pops_ht.Class)
        pop_ht = hl.read_table(ancestry_hybrid_ht_path(data_source, freeze))
        pop_joint_ht = pop_ht.join(regeneron_pops_ht, 'left')
        pop_joint_ht = pop_joint_ht.key_by(array_map=pop_joint_ht.s.split("_")[1])
        pop_joint_ht = pop_joint_ht.join(ukbb_ancestry, 'left')
        pop_joint_ht = pop_joint_ht.key_by('s')
        pop_joint_ht.write(get_joint_regeneron_ancestry_path(data_source, freeze), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

    parser.add_argument('--join_relatedness', help='Join relatedness table with Regeneron relatedness inference files.', action='store_true')
    parser.add_argument('--first_degree_kin_threshold', help='Upper and lower kinship threshold to use for first degree relatedness comma separated',
                        default="0.1767767,0.4")
    parser.add_argument('--second_degree_kin_threshold', help='Lower kinship threshold to use for second degree relatedness',
                        default=0.08838835)
    parser.add_argument('--ibd2_parent_offspring_threshold', help='IBD2 cutoff to determine parent offspring vs full sibling',
                        default=0.14)


    parser.add_argument('--filter_ukbb_ancestry', help='Subset UKBB phenotypes to self reported ancestry.', action='store_true')
    parser.add_argument('--join_ancestry', help='Join ancestry table with Regeneron and UKBB files.', action='store_true')


    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)