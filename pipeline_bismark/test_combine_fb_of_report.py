#! /usr/bin/python
import argparse
import os
from scipy.stats import binom_test
from statsmodels.sandbox.stats.multicomp import multipletests


def str2bool(v):
    # susendberg's function
    return v.lower() in ("yes", "true", "t", "1")


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


def get_rmet_type(crmet):
    if crmet <= 0.0:
        return 0  # 'unmethyed'
    elif 0.0 < crmet <= 0.30:
        return 1  # 'lowly'
    elif 0.30 < crmet <= 0.70:
        return 2  # 'moderately'
    elif 0.70 < crmet <= 1.00:
        return 3  # 'highly'


def combine_fb_of_deepreport(report_fp, cgposes):
    pos2info = {}
    for cgpos in cgposes:
        pos2info[cgpos] = [0.0, 0.0, 0, 0, 0, 0.0, 0, '']
    with open(report_fp, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split('\t')
            keytmp = (words[0], int(words[1]))
            if words[2] == '-':
                keytmp = (words[0], int(words[1]) - 1)
            prob0, prob1, met, unmet, coverage, cgtype = float(words[4]), float(words[5]), \
                int(words[6]), int(words[7]), int(words[8]), words[11]
            pos2info[keytmp][0] += prob0
            pos2info[keytmp][1] += prob1
            pos2info[keytmp][2] += met
            pos2info[keytmp][3] += unmet
            pos2info[keytmp][4] += coverage
            if words[2] == '+' or pos2info[keytmp][7] == '':
                pos2info[keytmp][7] = cgtype
    for cgpos in list(pos2info.keys()):
        if pos2info[cgpos][4] == 0:
            del pos2info[cgpos]
    mposinfo = []
    for cgpos in pos2info.keys():
        mposinfo.append(list(cgpos) + ['+', cgpos[1]] + pos2info[cgpos])
    mposinfo = sorted(mposinfo, key=lambda x: (x[0], x[1]))
    return mposinfo


def add_rmet_and_rtype(mposinfo):
    for mpos in mposinfo:
        mpos[9] = float(mpos[6]) / mpos[8]
        mpos[10] = get_rmet_type(mpos[9])
    return mposinfo


def combine_fb_of_bsreport(report_fp, cgposes):
    pos2info = {}
    for cgpos in cgposes:
        pos2info[cgpos] = [0, 0, 0, 0.0, 0.0, 1.0, 1.0, 0]
    with open(report_fp, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split('\t')
            keytmp = (words[0], int(words[1]))
            if words[2] == '-':
                keytmp = (words[0], int(words[1]) - 1)
            met, unmet, coverage = int(words[4]), int(words[5]), int(words[6])
            pos2info[keytmp][0] += met
            pos2info[keytmp][1] += unmet
            pos2info[keytmp][2] += coverage
    for cgpos in list(pos2info.keys()):
        if pos2info[cgpos][2] == 0:
            del pos2info[cgpos]
    mposinfo = []
    for cgpos in pos2info.keys():
        mposinfo.append(list(cgpos) + ['+', cgpos[1]] + pos2info[cgpos])
    mposinfo = sorted(mposinfo, key=lambda x: (x[0], x[1]))
    return mposinfo


def add_crmet_and_pqvalue(mposinfo, conversion_rate=0.9966, cal_pval=False, cal_qval=False):
    # cal rmet, cRmet, pvalue
    for mpos in mposinfo:
        mpos[7] = float(mpos[4]) / mpos[6]
        mpos[8] = (mpos[7] - (1 - conversion_rate)) / conversion_rate
        mpos[11] = get_rmet_type(mpos[8])
        if cal_pval:
            mpos[9] = binom_test(mpos[5], mpos[6], conversion_rate, alternative='less')
        else:
            mpos[9] = 0.0
    # cal qvalue
    if cal_qval:
        pvalues = [x[9] for x in mposinfo]
        qvalues = multipletests(pvalues, alpha=0.05, method='fdr_bh',
                                is_sorted=False, returnsorted=False)[1]
    else:
        qvalues = [x[9] for x in mposinfo]
    for i in range(0, len(mposinfo)):
        mposinfo[i][10] = qvalues[i]
    return mposinfo


def combine_fb_of_bsreport_encode(report_fp, cgposes):
    pos2info = {}
    for cgpos in cgposes:
        pos2info[cgpos] = [0, 0.0, 0.0]
    with open(report_fp, "r") as rf:
        # next(rf)
        for line in rf:
            words = line.strip().split('\t')
            keytmp = (words[0], int(words[1]))
            if words[5] == '-':
                keytmp = (words[0], int(words[1]) - 1)
            coverage, met = int(words[9]), float(words[10]) / 100 * int(words[9])
            try:
                pos2info[keytmp][0] += coverage
                pos2info[keytmp][1] += met
            except KeyError:
                pass
    for cgpos in list(pos2info.keys()):
        if pos2info[cgpos][0] == 0:
            del pos2info[cgpos]
        else:
            pos2info[cgpos][2] = pos2info[cgpos][1] / pos2info[cgpos][0]
    mposinfo = []
    for cgpos in pos2info.keys():
        mposinfo.append(list(cgpos) + pos2info[cgpos])
    mposinfo = sorted(mposinfo, key=lambda x: (x[0], x[1]))
    return mposinfo


def write_mpos2covinfo_bs(mclist, mcfile):

    with open(mcfile, 'w') as wf:
        wf.write('\t'.join(['chromosome', 'pos', 'strand', 'pos_in_strand',
                            'met', 'unmet', 'coverage', 'Rmet', 'cRmet',
                            'p.value', 'q.value', 'Rtype']) + '\n')
        for mctmp in mclist:
            wf.write('\t'.join(list(map(str, list(mctmp)))) + '\n')
    return mclist


def write_mpos2covinfo_deep(mclist, reportfp):
    with open(reportfp, 'w') as wf:
        wf.write('\t'.join(['chromosome', 'pos', 'strand', 'pos_in_strand', 'prob0', 'prob1',
                            'met', 'unmet', 'coverage', 'Rmet', 'rtype', 'cgtype']) + '\n')
        for mctmp in mclist:
            wf.write('\t'.join(list(map(str, list(mctmp)))) + '\n')
    return mclist


def write_mpos2covinfo_bs_encode(mclist, mcfile):
    with open(mcfile, 'w') as wf:
        wf.write('\t'.join(['chromosome', 'pos', 'coverage', 'met', 'Rmet']) + '\n')
        for mctmp in mclist:
            wf.write('\t'.join(list(map(str, list(mctmp)))) + '\n')
    return mclist


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--report_fp", help="the pred_in_ref report file path",
                        type=str, required=True)
    parser.add_argument('-t', "--rtype", help="bs or deepsignal, or bs_encode",
                        type=str, default='bs')
    parser.add_argument('--conversion_rate', type=float, default=1.0, required=False,
                        help='bisulfite conversion rate')
    parser.add_argument('-r', "--ref_fp", help="the directory of genome reference",
                        type=str, required=True)
    parser.add_argument('--contig', type=str, required=False, default='',
                        help='contig name need to be processed, if default, then all contigs are used')
    parser.add_argument("--cal_pval", type=str, default="no", help="cal qvalue or not")
    parser.add_argument("--cal_qval", type=str, default="no", help="cal qvalue or not")
    argv = parser.parse_args()

    report_fp = argv.report_fp
    rtype = argv.rtype
    ref_fp = argv.ref_fp
    contign = argv.contig
    conversion_rate = argv.conversion_rate
    cal_pval = str2bool(argv.cal_pval)
    cal_qval = str2bool(argv.cal_qval)

    print('start to get genome reference info..')
    refseq = DNAReference(ref_fp)
    contigname2contigseq = refseq.getcontigs()
    del refseq

    print('start to get motif poses in genome reference..')
    contig_cg_poses = set()
    if contign == '':
        for cgname in contigname2contigseq.keys():
            fcseq = contigname2contigseq[cgname]
            fposes = get_refloc_of_methysite_in_motif(fcseq, 'CG', 0)
            for fpos in fposes:
                contig_cg_poses.add((cgname, fpos))
    else:
        fcseq = contigname2contigseq[contign]
        fposes = get_refloc_of_methysite_in_motif(fcseq, 'CG', 0)
        for fpos in fposes:
            contig_cg_poses.add((contign, fpos))

    print('start to combine forward backward strands..')
    fname, fext = os.path.splitext(report_fp)
    wfp = fname + '.fb_combined' + fext
    if rtype == 'bs':
        mposinfo = combine_fb_of_bsreport(report_fp, contig_cg_poses)
        mposinfo = add_crmet_and_pqvalue(mposinfo, conversion_rate, cal_pval, cal_qval)
        write_mpos2covinfo_bs(mposinfo, wfp)
    elif rtype == 'deepsignal':
        mposinfo = combine_fb_of_deepreport(report_fp, contig_cg_poses)
        mposinfo = add_rmet_and_rtype(mposinfo)
        write_mpos2covinfo_deep(mposinfo, wfp)
    elif rtype == 'bs_encode':
        mposinfo = combine_fb_of_bsreport_encode(report_fp, contig_cg_poses)
        write_mpos2covinfo_bs_encode(mposinfo, wfp)
    else:
        pass


if __name__ == '__main__':
    main()
