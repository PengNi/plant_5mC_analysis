#! /usr/bin/python
import argparse
import os
import random
import numpy
import uuid

tmp_dir = "/home/nipeng"
abs_dir = os.path.dirname(os.path.realpath(__file__))
corr_rscript = "/".join([os.path.dirname(abs_dir), "/correlation_with_bs.cal_plot.general.R"])


def str2bool(v):
    # susendberg's function
    return v.lower() in ("yes", "true", "t", "1")


def read_r_stat(r_stat_fp, is_header=False):
    result = []
    with open(r_stat_fp, 'r') as rf:
        line = next(rf)
        if is_header:
            colnames = line.strip().split('\t')
        else:
            colnames = None
        words = next(rf).strip().split('\t')
        # todo: not safe to use -2
        for i in range(0, len(words)-2):
            result.append(int(words[i]))
        for i in range(len(words)-2, len(words)):
            result.append(float(words[i]))
    return result, colnames


def _read_one_mod_freq_file(freqfile):
    freqinfo = {}
    with open(freqfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            m_key = "\t".join([words[0], words[1], words[2]])
            pos_in_strand = words[3]
            methy_prob = float(words[4])
            unmethy_prob = float(words[5])
            methy_cov = int(words[6])
            unmethy_cov = int(words[7])
            cov = int(words[8])
            rmet = float(words[9])
            kmer = words[10]
            freqinfo[m_key] = [pos_in_strand, methy_prob, unmethy_prob, methy_cov, unmethy_cov,
                               cov, rmet, kmer]
    return freqinfo


def _get_combined_freq_file(freqfiles):
    freqinfo = {}
    freqkeys = set()
    for ffile in freqfiles:
        finfo_tmp = _read_one_mod_freq_file(ffile)
        for fkey in finfo_tmp.keys():
            if fkey not in freqkeys:
                freqkeys.add(fkey)
                freqinfo[fkey] = ["", 0.0, 0.0, 0, 0, 0, 0.0, ""]
            freqinfo[fkey][0] = finfo_tmp[fkey][0]
            freqinfo[fkey][1] += finfo_tmp[fkey][1]
            freqinfo[fkey][2] += finfo_tmp[fkey][2]
            freqinfo[fkey][3] += finfo_tmp[fkey][3]
            freqinfo[fkey][4] += finfo_tmp[fkey][4]
            freqinfo[fkey][5] += finfo_tmp[fkey][5]
            freqinfo[fkey][6] = freqinfo[fkey][3] / float(freqinfo[fkey][5])
            freqinfo[fkey][7] = finfo_tmp[fkey][7]
    return freqinfo


def _write_freqinfo(freqinfo, wfile):
    wf = open(wfile, "w")
    for fkey in freqinfo.keys():
        wstr = "\t".join([fkey, ] + list(map(str, freqinfo[fkey]))) + "\n"
        wf.write(wstr)
    wf.close()


def _run_rscript(freqfile, bsfile, bs_id):
    tmp_result_file = os.path.basename(freqfile) + ".stats"
    r_cmd = " ".join(['Rscript', corr_rscript, bsfile, freqfile, "bisulfite."+bs_id,
                      "deepsignal", "NCN", "no", tmp_dir.rstrip("/"), tmp_result_file])
    os.system(r_cmd)
    corr_stats, colnames = read_r_stat(tmp_dir + "/" + tmp_result_file, True)
    os.remove(tmp_dir + "/" + tmp_result_file)
    return corr_stats, colnames


def findsubsets(s, n):
    import itertools
    return list(itertools.combinations(s, n))


def eval_coverage_effect(args):
    wf = open(args.wfile, "w")

    modsfiles = args.modsfile
    modsfile_len = len(modsfiles)
    for modfile_idx in range(0, modsfile_len):
        wf.write("nano=={}, {}\n".format(modfile_idx, modsfiles[modfile_idx]))
    wf.write("\n")

    bsfiles = args.bsfile
    bsfile_len = len(bsfiles)
    for bs_idx in range(0, bsfile_len):
        wf.write("bs=={}, {}\n".format(bs_idx, bsfiles[bs_idx]))
    wf.write("\n")
    wf.flush()

    # coverage
    for coverage_idx in range(1, modsfile_len + 1):
        coverage_name = "{}0x".format(coverage_idx)
        combines = findsubsets([i for i in range(0, modsfile_len)], coverage_idx)
        random.shuffle(combines)
        iterations = args.repeat if args.repeat < len(combines) else len(combines)
        wf.write("================={}\n".format(coverage_name))
        # iteration
        stats_cove = []
        for iter_idx in range(0, iterations):
            # wf.write("=====iteration{}\n".format(iter_idx+1))
            iter_name = "_".join(list(map(str, list(combines[iter_idx]))))

            tmp_freqfile = tmp_dir + "/mod_freq_" + iter_name + "." + str(uuid.uuid1()) + ".tsv"
            tmp_modsfiles = [modsfiles[i] for i in combines[iter_idx]]
            tmp_freqinfo = _get_combined_freq_file(tmp_modsfiles)
            _write_freqinfo(tmp_freqinfo, tmp_freqfile)

            curid_corr_stats = []
            colnames = ""
            for bs_idx in range(0, len(bsfiles)):
                corr_stats, colnames = _run_rscript(tmp_freqfile, bsfiles[bs_idx], args.bs_id + str(bs_idx))
                curid_corr_stats.append(corr_stats)
            os.remove(tmp_freqfile)
            if iter_idx == 0:
                wf.write('\t'.join(['iter', ] + colnames) + '\n')
            for bs_idx in range(0, len(bsfiles)):
                curid_corr_tmp = curid_corr_stats[bs_idx]
                wf.write('\t'.join([str(iter_name) + "_vs_" + str(bs_idx), ] + list(map(str, curid_corr_tmp))) + '\n')
            curid_corr_stats = numpy.array(curid_corr_stats, dtype=numpy.float)
            curid_corr_stats_mean = numpy.mean(curid_corr_stats, 0)
            wf.write('\t'.join([str(iter_name) + "_bsmean", ] + list(map(str, curid_corr_stats_mean))) + '\n')
            stats_cove.append(curid_corr_stats_mean)
        stats_cove = numpy.array(stats_cove, dtype=numpy.float)
        stats_cove_mean = numpy.mean(stats_cove, 0)
        stats_cove_std = numpy.std(stats_cove, 0)
        wf.write('mean\t' + '\t'.join(list(map(str, stats_cove_mean))) + '\n')
        wf.write('std\t' + '\t'.join(list(map(str, stats_cove_std))) + '\n')
        wf.write("\n")
        wf.flush()
    wf.flush()
    wf.close()


def main():
    # python ~/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/performance_coverage_effect.py --modsfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.freq.tsv --modsfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.freq.tsv --modsfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.freq.tsv --modsfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.freq.tsv --modsfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.freq.tsv --bsfile bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --bsfile bs.poses/D1902827A-ZJ_364364.L1L2_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --bsfile bs.poses/D1902828A-ZJ_363363.L2L4_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --repeat 5 --wfile ~/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.diffcov_vs_bsrep123.tsv > ~/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.diffcov_vs_bsrep123.log 2>&1 &
    parser = argparse.ArgumentParser()
    parser.add_argument("--modsfile", action="append", type=str, required=True,
                        help="10x call_mods_freq file")
    parser.add_argument("--bsfile", type=str, action="append", required=True,
                        help="bs rmet file")
    parser.add_argument("--bs_id", type=str, required=False, default="rep",
                        help="replicate num, e.g.: rep1")
    parser.add_argument("--repeat", type=int, required=False, default=5,
                        help="random repeat times, default 5")
    parser.add_argument("--wfile", type=str, required=True,
                        help="")

    args = parser.parse_args()
    eval_coverage_effect(args)


if __name__ == '__main__':
    main()
