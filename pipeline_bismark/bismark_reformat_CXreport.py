#! /usr/bin/python
import os
import argparse
from scipy.stats import binom_test
from statsmodels.sandbox.stats.multicomp import multipletests


# J02 0.9967 J03 0.9966
methy_flags = {"Z", "X", "H", "U"}
unmethy_flags = {"z", "x", "h", "u"}
alphabeta = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
# colnames: 'chromosome', 'pos', 'strand', 'pos_in_strand', 'met', 'unmet', 'coverage',
# 'Rmet', 'cRmet', 'p.value', 'q.value', 'rtype'
# rtype: [0, 30], (30, 70], (70, 100]


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


def get_refloc_of_methysite_in_motif(seqstr, motif='CG', methyloc_in_motif=0):
    """

    :param seqstr:
    :param motif:
    :param methyloc_in_motif: 0-based
    :return:
    """
    strlen = len(seqstr)
    motiflen = len(motif)
    sites = []
    for i in range(0, strlen - motiflen + 1):
        if seqstr[i:i + motiflen] == motif:
            sites.append(i+methyloc_in_motif)
    return sites


def convert_motif_seq(ori_seq):
    alphabets = {'A': ['A', ], 'T': ['T', ], 'C': ['C', ], 'G': ['G', ],
                 'N': ['N', ], 'H': ['A', 'C', 'T']}
    outbases = []
    for bbase in ori_seq:
        outbases.append(alphabets[bbase])

    def recursive_permute(bases_list):
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


def get_all_motif_positions(contigs, contigname, motifs_str, meloc=0):
    motif_poses = []
    motifs = motifs_str.strip().split(',')

    motif_seqs = []
    for ori_motif in motifs:
        motif_seqs += convert_motif_seq(ori_motif)

    ctmpseq = contigs[contigname]
    ctmpseq_rc = complement_seq(ctmpseq)
    ctmplen = len(ctmpseq)

    # each motif
    poses = []
    for motif in motif_seqs:
        # forward strand
        f_poses = get_refloc_of_methysite_in_motif(ctmpseq, motif, meloc)
        poses += [(contigname, x, "+") for x in f_poses]
        # backward strand
        b_poses = get_refloc_of_methysite_in_motif(ctmpseq_rc, motif, meloc)
        b_poses = [ctmplen - 1 - x for x in b_poses]
        poses += [(contigname, x, '-') for x in b_poses]
    print('there are {} motifs in contig {}'.format(len(poses), contigname))
    if len(poses) > 0:
        motif_poses += sorted(poses, key=lambda p: (p[0], p[2], p[1]))
    return motif_poses


def get_rmet_type(crmet):
    if crmet <= 0.0:
        return 0  # 'unmethyed'
    elif 0.0 < crmet <= 0.30:
        return 1  # 'lowly'
    elif 0.30 < crmet <= 0.70:
        return 2  # 'moderately'
    elif 0.70 < crmet <= 1.00:
        return 3  # 'highly'


def get_mpositions_from_CXreport(cx_report_file, motifstr, conversion_rate=1.0, cal_pval=False):
    mpos2covinfo = {}
    with open(cx_report_file, 'r') as rf:
        for line in rf:
            words = line.strip().split('\t')
            keytmp = (words[0], int(words[1]), words[2])
            mcnt, unmcnt = int(words[3]), int(words[4])
            motif = words[5]
            if motif != motifstr:
                continue
            # if words[2] == "+":
            #     pos_strand = int(words[1])
            # else:
            #     pos_strand = len(contigs[words[0]]) - 1 - int(words[1])
            pos_strand = "-"
            if mcnt == 0 and unmcnt == 0:
                # ignore uncovered CGs
                # rtype = get_rmet_type(0.0)
                # p_value = binom_test(0, 0, conversion_rate, alternative='less')
                # mpos2covinfo[keytmp] = (pos_strand, mcnt, unmcnt, 0, 0, 0, p_value, 1, rtype)
                continue
            else:
                rmet = float(mcnt) / (mcnt + unmcnt)
                crmet = (rmet - (1 - conversion_rate)) / conversion_rate
                rtype = get_rmet_type(crmet)
                if cal_pval:
                    p_value = binom_test(unmcnt, mcnt + unmcnt, conversion_rate, alternative='less')
                else:
                    p_value = 0.0
                mpos2covinfo[keytmp] = (pos_strand, mcnt, unmcnt, mcnt + unmcnt,
                                        rmet, crmet, p_value, 1, rtype)
    return mpos2covinfo


def mposdict2list(mpos2covinfo):
    mclist = []
    for mpos in mpos2covinfo.keys():
        mclist.append(list(mpos) + list(mpos2covinfo[mpos]))
    mclist = sorted(mclist, key=lambda x: (x[0], x[1]))
    return mclist


def calculate_q_value(mclist):
    pvalues = [x[9] for x in mclist]
    try:
        qvalues = multipletests(pvalues, alpha=0.05, method='fdr_bh',
                                is_sorted=False, returnsorted=False)[1]
        return qvalues
    except ZeroDivisionError:
        return pvalues


def write_mpos2covinfo(mclist, mcfile):

    with open(mcfile, 'w') as wf:
        wf.write('\t'.join(['chromosome', 'pos', 'strand', 'pos_in_strand',
                            'met', 'unmet', 'coverage', 'Rmet', 'cRmet',
                            'p.value', 'q.value', 'Rtype']) + '\n')
        for mctmp in mclist:
            wf.write('\t'.join(list(map(str, list(mctmp)))) + '\n')
    return mclist


def write_mpos2covinfo_append(mclist, mcfile):

    with open(mcfile, 'a') as wf:
        for mctmp in mclist:
            wf.write('\t'.join(list(map(str, list(mctmp)))) + '\n')
    return mclist


def str2bool(v):
    # susendberg's function
    return v.lower() in ("yes", "true", "t", "1")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', "--motifs", default='CG',
                        help="motif",
                        type=str, required=False)
    parser.add_argument('--cx_report', type=str, default='', required=False,
                        help='CX_report file from Bismark')
    parser.add_argument('--conversion_rate', type=float, default=1.0, required=False,
                        help='bisulfite conversion rate')
    parser.add_argument("--cal_pval", type=str, default="no", help="cal qvalue or not")
    parser.add_argument("--cal_qval", type=str, default="no", help="cal qvalue or not")
    argv = parser.parse_args()
    motifstr = argv.motifs
    cxreport_fp = argv.cx_report
    conversion_rate = argv.conversion_rate

    fname, fext = os.path.splitext(cxreport_fp)
    motif_str = motifstr.strip().replace(',', '.')
    mcfp = fname + '.' + motif_str + fext
    wf = open(mcfp, 'w')
    wf.write('\t'.join(['chromosome', 'pos', 'strand', 'pos_in_strand',
                        'met', 'unmet', 'coverage', 'Rmet', 'cRmet',
                        'p.value', 'q.value', 'Rtype']) + '\n')
    wf.close()
    mpos2covinfo = get_mpositions_from_CXreport(cxreport_fp, motifstr, conversion_rate, str2bool(argv.cal_pval))
        
    mclist = mposdict2list(mpos2covinfo)
    if str2bool(argv.cal_qval):
        qvalues = calculate_q_value(mclist)
    else:
        qvalues = [x[9] for x in mclist]
    for i in range(0, len(mclist)):
        mclist[i][10] = qvalues[i]

    write_mpos2covinfo_append(mclist, mcfp)


if __name__ == '__main__':
    main()
