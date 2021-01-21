#! /usr/bin/python
import argparse
import os
from scipy import stats
import multiprocessing as mp
from multiprocessing import Queue
import time


def _cal_binom_pval(x, n, p=0.5):
    return stats.binom_test(x, n, p)


def fill_line_queue(input_path, queue_lines):
    print("read_process-{} started..".format(os.getpid()))
    with open(input_path, "r") as rf:
        next(rf)
        for line in rf:
            queue_lines.put(line)
            while queue_lines.qsize() > 100000:
                time.sleep(1)
    queue_lines.put("kill")
    print("read_process-{} finished..".format(os.getpid()))


#     1) chromorome
#     2) coordinate (1-based)
#     3) strand
#     4) sequence context (CG|CHG|CHH)
#     5) methylation ratio, calculated as #C_counts / #eff_CT_counts
#     6) number of effective total C+T counts on this locus (#eff_CT_counts)
#     7) number of total C counts on this locus (#C_counts)
#     8) number of total C+T counts on this locuso (#CT_counts)
#     9) number of total G counts on this locus of reverse strand (#rev_G_counts)
#     10) number of total G+A counts on this locus of reverse strand (#rev_GA_counts)
#     11) lower bound of 95% confidence interval of methylation ratio
#     12) upper bound of 95% confidence interval of methylation ratio
#     13) probability of observing the number of reads with either a C or a T by chance under the assumption
#         that the original distribution is 50%:50%
#     14) methylation call (5mC or C)
def _process_one_line(line):
    words = line.strip().split("\t")
    chrom, pos, strand, motif, rmet, effcov = words[0], int(words[1]) - 1, words[2], words[3], \
                                              float(words[4]), float(words[5])
    # prob = float(words[12])
    n_met = int(words[6])
    n_cov = int(round(effcov, 0))
    if n_met > n_cov:
        n_cov = n_met
    p_val = _cal_binom_pval(n_met, n_cov, 0.5)
    return motif, [chrom, str(pos), str(pos+1), str(p_val), str(effcov), strand, str(pos), str(pos + 1),
                   "0,0,0", str(effcov), str(int(round(rmet*100, 0)))]


def _process_line_queue(queue_lines, queue_write):
    print('cal_process-{} started..'.format(os.getpid()))
    while True:
        line = queue_lines.get()
        if line == "kill":
            queue_lines.put("kill")
            print('cal_process-{} finished..'.format(os.getpid()))
            break
        queue_write.put(_process_one_line(line))
        while queue_write.qsize() > 200000:
            time.sleep(1)


def _write_line(bsmapfile, queue_write):
    fname, fext = os.path.splitext(bsmapfile)

    wfile_cg = fname + ".methylbed." + "CG" + ".bed"
    wfile_chg = fname + ".methylbed." + "CHG" + ".bed"
    wfile_chh = fname + ".methylbed." + "CHH" + ".bed"
    wf_cg = open(wfile_cg, "w")
    wf_chg = open(wfile_chg, "w")
    wf_chh = open(wfile_chh, "w")
    print('write_process-{} started..'.format(os.getpid()))
    while True:
        # if queue_write.empty():
        #     time.sleep(1)
        #     continue
        item_write = queue_write.get()
        if item_write == "kill":
            print('write_process-{} finished'.format(os.getpid()))
            wf_cg.flush()
            wf_cg.close()
            wf_chg.flush()
            wf_chg.close()
            wf_chh.flush()
            wf_chh.close()
            break
        motif, strlist = item_write
        if motif == "CG":
            wf_cg.write("\t".join(strlist) + "\n")
        elif motif == "CHG":
            wf_chg.write("\t".join(strlist) + "\n")
        elif motif == "CHH":
            wf_chh.write("\t".join(strlist) + "\n")
        else:
            raise ValueError()


def convert_bsmap_file(args):

    queue_lines = Queue()
    queue_write = Queue()

    # read
    p_read = mp.Process(target=fill_line_queue, args=(args.bsmapfile, queue_lines))
    p_read.daemon = True
    p_read.start()

    p_cals = []
    nproc = args.nproc
    if nproc < 3:
        nproc = 3
    for _ in range(nproc - 2):
        p_cal = mp.Process(target=_process_line_queue, args=(queue_lines, queue_write))
        p_cal.daemon = True
        p_cal.start()
        p_cals.append(p_cal)

    p_w = mp.Process(target=_write_line, args=(args.bsmapfile, queue_write))
    p_w.daemon = True
    p_w.start()

    p_read.join()

    for p_cal in p_cals:
        p_cal.join()

    queue_write.put("kill")
    p_w.join()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bsmapfile", type=str, action="store", required=True,
                        help="bsmapfile")
    parser.add_argument("--nproc", "-p", action="store", type=int, default=3,
                        required=False,
                        help="number of processes to be used, default 3")
    args = parser.parse_args()

    convert_bsmap_file(args)


if __name__ == '__main__':
    main()
