#! /usr/bin/python
import argparse
from gff_reader import GFF3
from gff_reader import extract_region_by_attri
from gff_reader import extract_region_by_attri_conditionFeature
from gff_reader import get_kinds_of_a_attri
from gff_reader import get_kinds_of_a_attri_conditionFeature
from gff_reader import get_featurekinds_of_a_attrikey
from gff_reader import get_featurekinds_conditionIDstarts


def analyse_gff(args):
    gff3info = GFF3(args.gff3)
    gff3_eles = gff3info.get_eles()
    print("{} gff elements".format(len(gff3_eles)))

    features = gff3info.get_features()
    feature2counts = dict()
    for feature in features:
        feature2counts[feature] = 0
    for gff3_ele in gff3_eles:
        feature2counts[gff3_ele.get_feature()] += 1
    print("=====")
    for feature in sorted(feature2counts.keys()):
        print(feature, feature2counts[feature])
    print()

    print("=====")
    print("attri {} kinds:\n{}\n".format("locus_type",
                                         get_kinds_of_a_attri(gff3_eles,
                                                              "locus_type")))
    print("attri {} in feature-{} kinds:\n{}\n".format("locus_type", "gene",
                                                       get_kinds_of_a_attri_conditionFeature(gff3_eles,
                                                                                             "locus_type",
                                                                                             "gene")))
    print("attri {} in feature-{} kinds:\n{}\n".format("locus_type", "pseudogene",
                                                       get_kinds_of_a_attri_conditionFeature(gff3_eles,
                                                                                             "locus_type",
                                                                                             "pseudogene")))
    print("feature names of gffeles has attrikey-{}:\n{}".format("locus_type",
                                                                 get_featurekinds_of_a_attrikey(gff3_eles,
                                                                                                "locus_type")))

    print()
    print("{}=={}: {}".format("locus_type", "protein_coding",
                              len(extract_region_by_attri(gff3_eles, "locus_type", "protein_coding"))))
    print("{}=={}: {}".format("locus_type", "novel_transcribed_region",
                              len(extract_region_by_attri(gff3_eles, "locus_type", "novel_transcribed_region"))))
    print()
    print("{}=={}: {}".format("locus_type", "transposable_element_gene",
                              len(extract_region_by_attri(gff3_eles, "locus_type", "transposable_element_gene"))))
    print()
    print("{}=={}: {}".format("locus_type", "pseudogene",
                              len(extract_region_by_attri(gff3_eles, "locus_type", "pseudogene"))))

    print("\n=====")
    print("attri {} kinds:\n{}\n".format("biotype",
                                         get_kinds_of_a_attri(gff3_eles,
                                                              "biotype")))
    print("attri {} in feature-{} kinds:\n{}\n".format("biotype", "gene",
                                                       get_kinds_of_a_attri_conditionFeature(gff3_eles,
                                                                                             "biotype",
                                                                                             "gene")))
    print("attri {} in feature-{} kinds:\n{}\n".format("biotype", "pseudogene",
                                                       get_kinds_of_a_attri_conditionFeature(gff3_eles,
                                                                                             "biotype",
                                                                                             "pseudogene")))
    print("attri {} in feature-{} kinds:\n{}\n".format("biotype", "ncRNA_gene",
                                                       get_kinds_of_a_attri_conditionFeature(gff3_eles,
                                                                                             "biotype",
                                                                                             "ncRNA_gene")))

    print("feature names of gffeles has attrikey-{}:\n{}".format("biotype",
                                                                 get_featurekinds_of_a_attrikey(gff3_eles,
                                                                                                "biotype")))
    print("feature names of gffeles ID starts with-{}:\n{}".format("gene",
                                                                   get_featurekinds_conditionIDstarts(gff3_eles,
                                                                                                      "gene")))

    print()

    print("{}=={} in feature-{}: {}".format("biotype", "protein_coding", "gene",
                                            len(extract_region_by_attri_conditionFeature(gff3_eles, "biotype",
                                                                                         "protein_coding",
                                                                                         "gene"))))
    print("{}=={} in feature-{}: {}".format("biotype", "nontranslating_CDS", "gene",
                                            len(extract_region_by_attri_conditionFeature(gff3_eles, "biotype",
                                                                                         "nontranslating_CDS",
                                                                                         "gene"))))
    print()
    print("{}=={} in feature-{}: {}".format("biotype", "pseudogene", "pseudogene",
                                            len(extract_region_by_attri_conditionFeature(gff3_eles, "biotype",
                                                                                         "pseudogene",
                                                                                         "pseudogene"))))
    print("{}=={} in feature-{}: {}".format("biotype", "tRNA_pseudogene", "pseudogene",
                                            len(extract_region_by_attri_conditionFeature(gff3_eles, "biotype",
                                                                                         "tRNA_pseudogene",
                                                                                         "pseudogene"))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff3", type=str, required=True, help="")

    args = parser.parse_args()

    analyse_gff(args)


if __name__ == '__main__':
    main()
