#! /usr/bin/python
import argparse
from gff_reader import GFF3
from gff_reader import extract_region_by_feature
from gff_reader import extract_region_by_attri_conditionFeature
from gff_reader import extract_region_by_notattri_conditionFeature
import os

sep = "||"
# chromname_map_arab = {"NC_003070.9": "Chr1",
#                       "NC_003071.7": "Chr2",
#                       "NC_003074.8": "Chr3",
#                       "NC_003075.7": "Chr4",
#                       "NC_003076.8": "Chr5",
#                       "NC_037304.1": "ChrM",
#                       "NC_000932.1": "ChrC"}
chromname_map_arab = {"Chr1": "NC_003070.9",
                      "Chr2": "NC_003071.7",
                      "Chr3": "NC_003074.8",
                      "Chr4": "NC_003075.7",
                      "Chr5": "NC_003076.8",
                      "ChrM": "NC_037304.1",
                      "ChrC": "NC_000932.1"}

chroms_to_stat = {"arab": set(["NC_003070.9", "NC_003071.7", "NC_003074.8", "NC_003075.7", "NC_003076.8"]),
                  "rice": set(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"])}


# ======================
alphabeta = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def _alphabeta(letter):
    if letter in alphabeta.keys():
        return alphabeta[letter]
    return 'N'


def complement_seq(dnaseq):
    rdnaseq = dnaseq[::-1]
    comseq = ''
    try:
        comseq = ''.join([_alphabeta(x) for x in rdnaseq])
    except Exception:
        print('something wrong in the dna sequence.')
    return comseq


class DNAReference:
    def __init__(self, reffile):
        self._contignames = []
        self._contigs = {}  # contigname 2 contigseq
        with open(reffile, 'r') as rf:
            contigname = ''
            contigseq = ''
            for line in rf:
                if line.startswith('>'):
                    if contigname != '' and contigseq != '':
                        self._contigs[contigname] = contigseq
                        self._contignames.append(contigname)
                    contigname = line.strip()[1:].split(' ')[0]
                    contigseq = ''
                else:
                    # turn to upper case
                    contigseq += line.strip().upper()
            self._contigs[contigname] = contigseq
            self._contignames.append(contigname)

    def getcontigs(self):
        return self._contigs

    def getcontignames(self):
        return self._contignames


def get_refloc_of_methysite_in_motif(seqstr, motifset, methyloc_in_motif=0):
    """

    :param seqstr:
    :param motifset:
    :param methyloc_in_motif: 0-based
    :return:
    """
    motifset = set(motifset)
    strlen = len(seqstr)
    motiflen = len(list(motifset)[0])
    sites = []
    for i in range(0, strlen - motiflen + 1):
        if seqstr[i:i + motiflen] in motifset:
            sites.append(i+methyloc_in_motif)
    return sites


def convert_motif_seq(ori_seq):
    alphabets = {'A': ['A', ], 'T': ['T', ], 'C': ['C', ], 'G': ['G', ],
                 'N': ['N', ], 'H': ['A', 'C', 'T']}
    outbases = []
    for bbase in ori_seq:
        outbases.append(alphabets[bbase])

    def recursive_permute(bases_list):
        if len(bases_list) == 1:
            return bases_list[0]
        if len(bases_list) == 2:
            pseqs = []
            for fbase in bases_list[0]:
                for sbase in bases_list[1]:
                    pseqs.append(fbase + sbase)
            return pseqs
        else:
            pseqs = recursive_permute(bases_list[1:])
            pseq_list = [bases_list[0], pseqs]
            return recursive_permute(pseq_list)
    return recursive_permute(outbases)


def get_motif_seqs(motifs):
    ori_motif_seqs = motifs.strip().split(',')

    motif_seqs = []
    for ori_motif in ori_motif_seqs:
        motif_seqs += convert_motif_seq(ori_motif.strip().upper())
    return motif_seqs


# read gff3===================================
def extract_bed6_from_geles(geles):
    regions = []
    for gele in geles:
        chrom, strand, start, end = gele.get_chromosome(), gele.get_strand(), gele.get_start(), gele.get_end()
        gid = gele.get_id()
        regions.append([chrom, start, end, gid, 0, strand])
    return regions


def change_arab_chromname(regions):
    regions_new = [[chromname_map_arab[region[0]], ] + region[1:] for region in regions]
    return regions_new


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
    regions_pcg = change_arab_chromname(extract_bed6_from_geles(geles_pcg))
    regions_pcg = [x for x in regions_pcg if x[0] in chroms_to_stat["arab"]]

    regions_ncg = change_arab_chromname(extract_bed6_from_geles(geles_ncg))
    regions_ncg = [x for x in regions_ncg if x[0] in chroms_to_stat["arab"]]

    regions_pg = change_arab_chromname(extract_bed6_from_geles(geles_pg))
    regions_pg = [x for x in regions_pg if x[0] in chroms_to_stat["arab"]]

    regions_teg = change_arab_chromname(extract_bed6_from_geles(geles_teg))
    regions_teg = [x for x in regions_teg if x[0] in chroms_to_stat["arab"]]

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
    regions_pcg = extract_bed6_from_geles(geles_pcg)
    regions_pcg = [x for x in regions_pcg if x[0] in chroms_to_stat["rice"]]

    regions_ncg = extract_bed6_from_geles(geles_ncg)
    regions_ncg = [x for x in regions_ncg if x[0] in chroms_to_stat["rice"]]

    regions_pg = extract_bed6_from_geles(geles_pg)
    regions_pg = [x for x in regions_pg if x[0] in chroms_to_stat["rice"]]

    return {"protein-coding genes": regions_pcg, "non-coding genes": regions_ncg,
            "pseudogenes": regions_pg}


def generate_transcript_regions_from_gff3(gff3_eles, species):
    geles_5utr = extract_region_by_feature(gff3_eles, "five_prime_UTR")
    geles_cds = extract_region_by_feature(gff3_eles, "CDS")
    geles_3utr = extract_region_by_feature(gff3_eles, "three_prime_UTR")
    print("==5'UTR: {}; CDS: {}; 3'UTR: {}".format(len(geles_5utr),
                                                   len(geles_cds),
                                                   len(geles_3utr)))
    regions_5utr = extract_bed6_from_geles(geles_5utr)
    if species == "arab":
        regions_5utr = change_arab_chromname(regions_5utr)
    regions_5utr = [x for x in regions_5utr if x[0] in chroms_to_stat[species]]

    regions_cds = extract_bed6_from_geles(geles_cds)
    if species == "arab":
        regions_cds = change_arab_chromname(regions_cds)
    regions_cds = [x for x in regions_cds if x[0] in chroms_to_stat[species]]

    regions_3utr = extract_bed6_from_geles(geles_3utr)
    if species == "arab":
        regions_3utr = change_arab_chromname(regions_3utr)
    regions_3utr = [x for x in regions_3utr if x[0] in chroms_to_stat[species]]

    return {"5'UTR": regions_5utr, "CDS": regions_cds,
            "3'UTR": regions_3utr}


def read_cytosine_cov_info_from_rmetfile(rmetfile):
    sites = set()
    sites2cov = dict()
    with open(rmetfile, "r") as rf:
        if ".CX_report." in rmetfile:
            next(rf)
        for line in rf:
            words = line.strip().split("\t")
            if rmetfile.endswith(".bed"):
                chrom, pos, strand, cov = words[0], int(words[1]), words[5], int(words[9])
            elif rmetfile.endswith(".freq.tsv"):
                chrom, pos, strand, cov = words[0], int(words[1]), words[2], int(words[8])
            elif ".CX_report." in rmetfile:
                chrom, pos, strand, cov = words[0], int(words[1]), words[2], int(words[6])
            else:
                raise ValueError()
            sites.add(sep.join([chrom, str(pos), strand]))
            sites2cov[sep.join([chrom, str(pos), strand])] = cov
    return sites, sites2cov


def gene_cover_stats(gene_regions, sites, sites2cov, contigs, motifstr, mod_loc, cov_cf):
    motifs = get_motif_seqs(motifstr)
    regionsinfo = []

    for genetype in gene_regions.keys():
        gene_region = gene_regions[genetype]

        for region in gene_region:
            chrom, start, end, name, score, strand = region
            if strand == "+":
                rep_cs = get_refloc_of_methysite_in_motif(contigs[chrom][start:end], motifs, mod_loc)
                motifsites = [x + start for x in rep_cs]
            else:
                rep_cs = get_refloc_of_methysite_in_motif(complement_seq(contigs[chrom][start:end]),
                                                          motifs, mod_loc)
                motifsites = [end - 1 - x for x in rep_cs]
            if len(motifsites) == 0:
                continue
            n_cov1, n_covcf = 0, 0
            for msite in motifsites:
                key = sep.join([chrom, str(msite), strand])
                if key in sites:
                    coverage = sites2cov[key]
                    if coverage >= 1:
                        n_cov1 += 1
                    if coverage >= cov_cf:
                        n_covcf += 1
            regionsinfo.append([chrom, start, end, name, score, strand, end-start,
                                motifstr, len(motifsites), n_covcf,
                                n_covcf/float(len(motifsites)), n_cov1, n_cov1/float(len(motifsites)),
                                genetype])
    return regionsinfo


def main():
    parser = argparse.ArgumentParser("extract gene region locs of 4 classes: protein-coding; "
                                     "non-coding; transposable_element_gene (transposon); pseudogene")
    parser.add_argument("--gff3", type=str, required=True, help="")
    parser.add_argument("--cytosine_info", type=str, required=True,
                        action="append",
                        help="cytosine rmet (nano/bs) file, best 50x+,"
                             "order as motifs")
    parser.add_argument("--cov_cf", type=int, default=5, required=False,
                        help="")
    parser.add_argument("--technique", type=str, default="",
                        required=False)

    parser.add_argument("--ref", type=str, required=True,
                        help="")
    parser.add_argument('-m', "--motifs", default='CG',
                        help="motifs, splited by comma",
                        type=str, required=False)
    parser.add_argument('--mloc_in_motif', type=int, required=False,
                        default=0,
                        help='0-based location of the methylation base in the motif, default 0')

    parser.add_argument("--species", type=str, required=True, choices=["arab", "rice"], help="")

    parser.add_argument("--wfile", type=str, required=True, default=None)

    args = parser.parse_args()

    technique = args.technique
    dnaref = DNAReference(args.ref)
    dnacontigs = dnaref.getcontigs()

    gff3info = GFF3(args.gff3)
    gff3_eles = gff3info.get_eles()
    print("{} gff elements from file-{}".format(len(gff3_eles), args.gff3))

    if args.species == "arab":
        genes = generate_gene_regions_from_gff3_arab(gff3_eles)
    elif args.species == "rice":
        genes = generate_gene_regions_from_gff3_rice(gff3_eles)
    else:
        raise ValueError()
    trans = generate_transcript_regions_from_gff3(gff3_eles, args.species)

    open(args.wfile, "w").close()
    motifs = args.motifs.split(",")
    cytosine_info_files = args.cytosine_info
    for i in range(len(cytosine_info_files)):
        cytosine_info_file = cytosine_info_files[i]
        motifstr = motifs[i]

        csites, csites2cov = read_cytosine_cov_info_from_rmetfile(cytosine_info_file)

        genecovinfo = gene_cover_stats(genes, csites, csites2cov, dnacontigs,
                                       motifstr, args.mloc_in_motif, args.cov_cf)
        trancovinfo = gene_cover_stats(trans, csites, csites2cov, dnacontigs,
                                       motifstr, args.mloc_in_motif, args.cov_cf)

        wf = open(args.wfile, "a")
        for ginfo in [genecovinfo, trancovinfo]:
            for gitem in ginfo:
                wf.write("\t".join(list(map(str, gitem))) + "\t" + technique + "\n")
        wf.close()


if __name__ == '__main__':
    main()
