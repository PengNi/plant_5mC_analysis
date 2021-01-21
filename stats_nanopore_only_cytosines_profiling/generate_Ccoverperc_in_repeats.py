#! /usr/bin/python
import argparse
import os

sep = "||"
chromname_map_arab = {"NC_003070.9": "Chr1",
                      "NC_003071.7": "Chr2",
                      "NC_003074.8": "Chr3",
                      "NC_003075.7": "Chr4",
                      "NC_003076.8": "Chr5",
                      "NC_037304.1": "ChrM",
                      "NC_000932.1": "ChrC"}
chromname_map_rice = {"NC_029256.1": "1",
                      "NC_029257.1": "2",
                      "NC_029258.1": "3",
                      "NC_029259.1": "4",
                      "NC_029260.1": "5",
                      "NC_029261.1": "6",
                      "NC_029262.1": "7",
                      "NC_029263.1": "8",
                      "NC_029264.1": "9",
                      "NC_029265.1": "10",
                      "NC_029266.1": "11",
                      "NC_029267.1": "12"}

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


# ====analysis
def stat_posinfo(posinfo):
    motif2num = {}
    for tpos in posinfo:
        chrom, loc, strand, motif = tpos
        if motif not in motif2num.keys():
            motif2num[motif] = 0
        motif2num[motif] += 1
    return motif2num


# read trf .dat file (1-based, [a, b])
# TODO: maybe not right, but now, we take the following 3 regions as 3 indepdent tandem repeats:
# TODO: inverted repeats are precessed the same as tandem repeats as well.
# 640415 640446 16 2.0 16 93 0 55 43 9 6 40 1.62
# 640240 641187 185 4.6 203 81 12 1177 44 8 13 33 1.75
# 639444 641187 387 4.4 389 84 10 1963 44 8 12 34 1.74
def read_regioninfo_from_trfdatfile(datfile, species="arab"):
    regions = []
    with open(datfile, 'r') as rf:
        chromname = ""
        count = 0
        for line in rf:
            if line.startswith("Sequence:"):
                chromname = line.strip().split(" ")[1]
            elif line[0].isdigit():
                words = line.strip().split(" ")
                start, end = int(words[0]) - 1, int(words[1])
                regions.append((chromname, start, end, "tr"+str(count), "0", "+"))
                regions.append((chromname, start, end, "tr"+str(count), "0", "-"))
                count += 1
    regions = [x for x in regions if x[0] in chroms_to_stat[species]]
    return regions


# read irf .dat file (1-based, [a, b])
def read_regioninfo_from_irfdatfile(datfile, species="arab"):
    regions = []
    with open(datfile, 'r') as rf:
        chromname = ""
        count = 0
        for line in rf:
            if line.startswith("Sequence:"):
                chromname = line.strip().split(" ")[1]
            elif line[0].isdigit():
                words = line.strip().split(" ")
                start1, end1 = int(words[0]) - 1, int(words[1])
                start2, end2 = int(words[3]) - 1, int(words[4])
                # regions.append((chromname, start1, end1, "ir"+str(count)+"_l", "0", "+"))
                # regions.append((chromname, start1, end1, "ir"+str(count)+"_l", "0", "-"))
                # regions.append((chromname, start2, end2, "ir"+str(count)+"_r", "0", "+"))
                # regions.append((chromname, start2, end2, "ir"+str(count)+"_r", "0", "-"))
                regions.append(((chromname, start1, end1, "ir"+str(count)+"_l", "0", "+"),
                                (chromname, start2, end2, "ir" + str(count) + "_r", "0", "+")))
                regions.append(((chromname, start1, end1, "ir"+str(count)+"_l", "0", "-"),
                                (chromname, start2, end2, "ir" + str(count) + "_r", "0", "-")))
                count += 1
    regions = [x for x in regions if x[0][0] in chroms_to_stat[species]]
    return regions


def read_regioninfo_from_bedfile(bedfile, species):
    regions = []
    with open(bedfile, 'r') as rf:
        for line in rf:
            if not line.startswith("#"):
                words = line.strip().split("\t")
                chrom, start, end, name, score, strand = words[0], int(words[1]), int(words[2]), words[3], \
                    int(words[4]), words[5]
                if species == "rice":
                    try:
                        chrom = chromname_map_rice[chrom]
                    except KeyError:
                        chrom = chrom
                regions.append((chrom, start, end, name, score, strand))
    regions = [x for x in regions if x[0] in chroms_to_stat[species]]
    return regions


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


def get_numinfo_of_repeat(repeat, sites, sites2cov, contigs, motifs, mod_loc, cov_cf):
    chrom, start, end, name, score, strand = repeat
    if strand == ".":
        strand = "+"
    if strand == "+":
        rep_cs = get_refloc_of_methysite_in_motif(contigs[chrom][start:end], motifs, mod_loc)
        motifsites = [x + start for x in rep_cs]
    else:
        rep_cs = get_refloc_of_methysite_in_motif(complement_seq(contigs[chrom][start:end]),
                                                  motifs, mod_loc)
        motifsites = [end - 1 - x for x in rep_cs]
    n_cov1, n_covcf = 0, 0
    if len(motifsites) == 0:
        return chrom, start, end, name, score, strand, 0, 0, 0
    for msite in motifsites:
        key = sep.join([chrom, str(msite), strand])
        if key in sites:
            coverage = sites2cov[key]
            if coverage >= 1:
                n_cov1 += 1
            if coverage >= cov_cf:
                n_covcf += 1
    return chrom, start, end, name, score, strand, len(motifsites), n_covcf, n_cov1


def repeat_cover_stats(repeats, sites, sites2cov, contigs, motifstr, mod_loc, cov_cf, repeattype):
    motifs = get_motif_seqs(motifstr)
    repeatsinfo = []
    for repeat in repeats:
        if len(repeat) == 2:  # irf
            repeat1, repeat2 = repeat
            repeat_covinfo1 = get_numinfo_of_repeat(repeat1, sites, sites2cov, contigs, motifs, mod_loc, cov_cf)
            repeat_covinfo2 = get_numinfo_of_repeat(repeat2, sites, sites2cov, contigs, motifs, mod_loc, cov_cf)
            motifnum, n_covcf, n_cov1, repeat_len = 0, 0, 0, 0
            chrom1, start1, end1, name1, score1, strand1, motifnum1, n_covcf1, n_cov11 = repeat_covinfo1
            chrom2, start2, end2, name2, score2, strand2, motifnum2, n_covcf2, n_cov12 = repeat_covinfo2

            motifnum += motifnum1
            n_covcf += n_covcf1
            n_cov1 += n_cov11
            repeat_len += end1 - start1

            motifnum += motifnum2
            n_covcf += n_covcf2
            n_cov1 += n_cov12
            repeat_len += end2 - start2

            if motifnum > 0:
                repeatsinfo.append([chrom1, start1, end2, name1, score1, strand1, repeat_len, motifstr,
                                    motifnum, n_covcf, n_covcf/float(motifnum), n_cov1,
                                    n_cov1 / float(motifnum), repeattype])
        else:
            repeat_covinfo = get_numinfo_of_repeat(repeat, sites, sites2cov, contigs, motifs, mod_loc, cov_cf)
            chrom, start, end, name, score, strand, motifnum, n_covcf, n_cov1 = repeat_covinfo
            if motifnum > 0:
                repeat_len = end - start
                repeatsinfo.append([chrom, start, end, name, score, strand, repeat_len, motifstr,
                                    motifnum, n_covcf, n_covcf/float(motifnum), n_cov1,
                                    n_cov1/float(motifnum), repeattype])
    return repeatsinfo


def stat_cytosine_ratio_in_repeats(args):
    technique = args.technique
    dnaref = DNAReference(args.ref)
    dnacontigs = dnaref.getcontigs()

    repeat_window = read_regioninfo_from_bedfile(args.windowmasker_bed, args.species)
    # print("windowmasker region num: {}".format(len(repeat_window)))
    repeat_repeat = read_regioninfo_from_bedfile(args.repeatmasker_bed, args.species)
    # print("repeatmasker region num: {}".format(len(repeat_repeat)))
    repeat_trf = read_regioninfo_from_trfdatfile(args.trf_dat, args.species)
    repeat_irf = read_regioninfo_from_irfdatfile(args.irf_dat, args.species)

    open(args.wfile, "w").close()
    motifs = args.motifs.split(",")
    cytosine_info_files = args.cytosine_info
    for i in range(len(cytosine_info_files)):
        cytosine_info_file = cytosine_info_files[i]
        motifstr = motifs[i]

        csites, csites2cov = read_cytosine_cov_info_from_rmetfile(cytosine_info_file)

        repeatinfo_window = repeat_cover_stats(repeat_window, csites, csites2cov, dnacontigs,
                                               motifstr, args.mloc_in_motif, args.cov_cf,
                                               "Repeats(WindowMasker)")
        print("Repeats(WindowMasker) len: {}".format(len(repeatinfo_window)))
        repeatinfo_repeat = repeat_cover_stats(repeat_repeat, csites, csites2cov, dnacontigs,
                                               motifstr, args.mloc_in_motif, args.cov_cf,
                                               "Repeats(RepeatMasker)")
        print("Repeats(RepeatMasker) len: {}".format(len(repeatinfo_repeat)))
        repeatinfo_trf = repeat_cover_stats(repeat_trf, csites, csites2cov, dnacontigs,
                                            motifstr, args.mloc_in_motif, args.cov_cf,
                                            "Tandem repeats")
        print("Tandem repeats len: {}".format(len(repeatinfo_trf)))
        repeatinfo_irf = repeat_cover_stats(repeat_irf, csites, csites2cov, dnacontigs,
                                            motifstr, args.mloc_in_motif, args.cov_cf,
                                            "Inverted repeats")
        print("Inverted repeats len: {}".format(len(repeatinfo_irf)))

        wf = open(args.wfile, "a")
        for repeatinfo in [repeatinfo_window, repeatinfo_repeat, repeatinfo_trf, repeatinfo_irf]:
            for repeatitem in repeatinfo:
                wf.write("\t".join(list(map(str, repeatitem))) + "\t" + technique + "\n")
        wf.close()


def main():
    parser = argparse.ArgumentParser("repeat regions")
    parser.add_argument("--cytosine_info", type=str, required=True,
                        action="append",
                        help="cytosine rmet (nano/bs) file, best 50x+, order as motifs")
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

    parser.add_argument("--species", type=str, required=False, default="arab", choices=["arab", "rice"], help="")
    parser.add_argument("--windowmasker_bed", type=str, required=True,
                        help="")
    parser.add_argument("--repeatmasker_bed", type=str, required=True,
                        help="")
    parser.add_argument("--trf_dat", type=str, required=True,
                        help="")
    parser.add_argument("--irf_dat", type=str, required=True,
                        help="")

    parser.add_argument("--wfile", type=str, required=True, default=None)

    args = parser.parse_args()
    stat_cytosine_ratio_in_repeats(args)


if __name__ == '__main__':
    main()
