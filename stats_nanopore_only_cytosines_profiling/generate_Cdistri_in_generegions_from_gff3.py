#! /usr/bin/python
import argparse
from gff_reader import GFF3
from gff_reader import extract_region_by_feature
from gff_reader import extract_region_by_attri_conditionFeature
from gff_reader import extract_region_by_notattri_conditionFeature
import os

sep = "||"
chromname_map_arab = {"NC_003070.9": "Chr1",
                      "NC_003071.7": "Chr2",
                      "NC_003074.8": "Chr3",
                      "NC_003075.7": "Chr4",
                      "NC_003076.8": "Chr5",
                      "NC_037304.1": "ChrM",
                      "NC_000932.1": "ChrC"}


def _combine_regions(gff3_eles):
    chrom2regions = dict()
    for gele in gff3_eles:
        chrom, strand, start, end = gele.get_chromosome(), gele.get_strand(), gele.get_start(), gele.get_end()
        keytmp = sep.join([chrom, strand])
        if keytmp not in chrom2regions.keys():
            chrom2regions[keytmp] = []
        chrom2regions[keytmp].append((start, end))
    return chrom2regions


def generate_gene_regions_from_gff3_arab(gff3_eles):
    # Araport11_GFF3_genes_transposons.201606.gff
    # protein-coding genes; non-coding genes; pseudogenes; transposons (transposable_element_genes)
    geles_pcg = extract_region_by_attri_conditionFeature(gff3_eles, "locus_type", "protein_coding", "gene")
    geles_ncg = extract_region_by_notattri_conditionFeature(gff3_eles, "locus_type", "protein_coding", "gene")
    geles_pg = extract_region_by_attri_conditionFeature(gff3_eles, "locus_type", "pseudogene", "pseudogene")
    geles_teg = extract_region_by_attri_conditionFeature(gff3_eles, "locus_type", "transposable_element_gene",
                                                         "transposable_element_gene")
    print("==protein-coding: {}; non-coding: {}; pseudogene: {}; transposons: {}".format(len(geles_pcg),
                                                                                         len(geles_ncg),
                                                                                         len(geles_pg),
                                                                                         len(geles_teg)))
    regions_pcg = _combine_regions(geles_pcg)
    regions_ncg = _combine_regions(geles_ncg)
    regions_pg = _combine_regions(geles_pg)
    regions_teg = _combine_regions(geles_teg)
    return {"protein-coding genes": regions_pcg, "non-coding genes": regions_ncg,
            "pseudogenes": regions_pg, "transposons": regions_teg}


def generate_gene_regions_from_gff3_rice(gff3_eles):
    # Oryza_sativa.IRGSP-1.0.45.gff3
    # protein-coding genes; non-coding genes; pseudogenes
    geles_pcg = extract_region_by_feature(gff3_eles, "gene")
    geles_ncg = extract_region_by_feature(gff3_eles, "ncRNA_gene")
    geles_pg = extract_region_by_feature(gff3_eles, "pseudogene")
    print("==protein-coding: {}; non-coding: {}; pseudogene: {}".format(len(geles_pcg),
                                                                        len(geles_ncg),
                                                                        len(geles_pg)))
    regions_pcg = _combine_regions(geles_pcg)
    regions_ncg = _combine_regions(geles_ncg)
    regions_pg = _combine_regions(geles_pg)
    return {"protein-coding genes": regions_pcg, "non-coding genes": regions_ncg,
            "pseudogenes": regions_pg}


def generate_transcript_regions_from_gff3(gff3_eles):
    geles_5utr = extract_region_by_feature(gff3_eles, "five_prime_UTR")
    geles_cds = extract_region_by_feature(gff3_eles, "CDS")
    geles_3utr = extract_region_by_feature(gff3_eles, "three_prime_UTR")
    print("==5'UTR: {}; CDS: {}; 3'UTR: {}".format(len(geles_5utr),
                                                   len(geles_cds),
                                                   len(geles_3utr)))
    regions_5utr = _combine_regions(geles_5utr)
    regions_cds = _combine_regions(geles_cds)
    regions_3utr = _combine_regions(geles_3utr)
    return {"5'UTR": regions_5utr, "CDS": regions_cds,
            "3'UTR": regions_3utr}


def read_cytosine_info(cytosine_info_file):
    posinfo = []
    with open(cytosine_info_file, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            chrom, loc, strand, motif = words[0], int(words[1]), words[2], words[-1]
            posinfo.append((chrom, loc, strand, motif))
    return posinfo


def stat_posinfo(posinfo):
    motif2num = {}
    for tpos in posinfo:
        chrom, loc, strand, motif = tpos
        if motif not in motif2num.keys():
            motif2num[motif] = 0
        motif2num[motif] += 1
    return motif2num


def reformat_posinfo_for_cmp(posinfo, species):
    chrom2motifsites = {}
    for tposinfo in posinfo:
        chrom, loc, strand, motif = tposinfo
        if species == "arab":
            keytmp = sep.join([chromname_map_arab[chrom], strand])
        else:
            keytmp = sep.join([chrom, strand])
        if keytmp not in chrom2motifsites.keys():
            chrom2motifsites[keytmp] = {}
        if motif not in chrom2motifsites[keytmp].keys():
            chrom2motifsites[keytmp][motif] = []
        chrom2motifsites[keytmp][motif].append(loc)
    return chrom2motifsites


def get_pos_distribution_in_generegion(rgtype2regions, chrom2motifsites):
    rgtype2stats = {}
    for rgtype in rgtype2regions.keys():
        motif2count = {}
        chrom2regions = rgtype2regions[rgtype]
        for chromkey in chrom2regions.keys():
            if chromkey in chrom2motifsites.keys():
                # print(rgtype, chromkey)
                regionsites = set()
                for region in chrom2regions[chromkey]:
                    start, end = region
                    for locidx in range(start, end):
                        regionsites.add(locidx)

                motif2sites = chrom2motifsites[chromkey]
                for motif in motif2sites.keys():
                    sitestmp = set(motif2sites[motif])
                    num_sites_in_region = len(sitestmp.intersection(regionsites))
                    if motif not in motif2count.keys():
                        motif2count[motif] = 0
                    motif2count[motif] += num_sites_in_region
        rgtype2stats[rgtype] = motif2count
    return rgtype2stats


def main():
    # python generate_Cdistri_in_generegions_from_gff3.py --gff3 Araport11_GFF3_genes_transposons.201606.sorted.gff
    # --cytosine_info ninanjie.c_count.bs_3reps_vs_nano50x.only_nanopore_detected_Cs.pos.txt --species arab
    parser = argparse.ArgumentParser("extract gene region locs of 4 classes: protein-coding; non-coding; "
                                     "transposable_element_gene (transposon); pseudogene")
    parser.add_argument("--gff3", type=str, required=True, help="")
    parser.add_argument("--cytosine_info", type=str, required=True,
                        help="cytosine location file, from generate_nanopore_only_cytosines_info.py")
    parser.add_argument("--species", type=str, required=True, choices=["arab", "rice"], help="")

    args = parser.parse_args()

    gff3info = GFF3(args.gff3)
    gff3_eles = gff3info.get_eles()
    print("{} gff elements from file-{}".format(len(gff3_eles), args.gff3))

    fname, fext = os.path.splitext(args.cytosine_info)
    posinfo = read_cytosine_info(args.cytosine_info)
    motif2counts = stat_posinfo(posinfo)
    print(motif2counts)
    cytosines_nano = reformat_posinfo_for_cmp(posinfo, args.species)
    totalnum = len(posinfo)

    # stat1 ====================================
    if args.species == "arab":
        rgtype2regions = generate_gene_regions_from_gff3_arab(gff3_eles)
    elif args.species == "rice":
        rgtype2regions = generate_gene_regions_from_gff3_rice(gff3_eles)
    else:
        raise ValueError("--species is not set rightly")

    rgtype2stats = get_pos_distribution_in_generegion(rgtype2regions, cytosines_nano)
    wfile1 = fname + ".cnum_in_gene.txt"
    with open(wfile1, "w") as wf:
        wf.write("\t".join(["region", "motif", "moitf_num", "motif_total", "c_total"]) + "\n")
        for rgtype in rgtype2stats.keys():
            for motif in rgtype2stats[rgtype].keys():
                wf.write("\t".join([rgtype, motif, str(rgtype2stats[rgtype][motif]),
                                    str(motif2counts[motif]), str(totalnum)]) + "\n")

    # stat2 ====================================
    rgtype2regions = generate_transcript_regions_from_gff3(gff3_eles)
    rgtype2stats = get_pos_distribution_in_generegion(rgtype2regions, cytosines_nano)
    wfile2 = fname + ".cnum_in_transcript.txt"
    with open(wfile2, "w") as wf:
        wf.write("\t".join(["region", "motif", "moitf_num", "motif_total", "c_total"]) + "\n")
        for rgtype in rgtype2stats.keys():
            for motif in rgtype2stats[rgtype].keys():
                wf.write("\t".join([rgtype, motif, str(rgtype2stats[rgtype][motif]),
                                    str(motif2counts[motif]), str(totalnum)]) + "\n")


if __name__ == '__main__':
    main()
