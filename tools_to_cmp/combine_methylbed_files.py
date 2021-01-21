#! /usr/bin/python
import argparse
import os
import random
import numpy
import uuid
import pandas as pd
import scipy.stats
import numpy as np


def str2bool(v):
    # susendberg's function
    return v.lower() in ("yes", "true", "t", "1")


def _read_one_mod_freq_file(freqfile):
    # methylbed format
    # "chromosome", "pos", "end", "na1", "na2", "strand", "na3", "na4", "na5", "coverage", "rpercent"
    freqinfo = {}
    with open(freqfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            m_key = "\t".join([words[0], words[1], words[5]])
            cov = float(words[9])
            rmet = float(words[10]) / 100
            methy_cov = rmet * cov
            freqinfo[m_key] = [methy_cov, cov, rmet]
    return freqinfo


def _get_combined_freq_info(bedfiles):
    freqinfo = {}
    freqkeys = set()
    for bfile in bedfiles:
        finfo = _read_one_mod_freq_file(bfile)
        for fkey in finfo.keys():
            if fkey not in freqkeys:
                freqkeys.add(fkey)
                freqinfo[fkey] = [0.0, 0.0, 0.0]
            freqinfo[fkey][0] += finfo[fkey][0]
            freqinfo[fkey][1] += finfo[fkey][1]
            if freqinfo[fkey][1] > 0:
                freqinfo[fkey][2] = freqinfo[fkey][0] / float(freqinfo[fkey][1])
    return freqinfo


def _write_freqinfo(freqinfo, wfile):
    wf = open(wfile, "w")
    for fkey in freqinfo.keys():
        key_items = fkey.split("\t")
        # "chromosome", "pos", "end", "na1", "na2", "strand", "na3", "na4", "na5", "coverage", "rpercent"
        chrom, pos, strand = key_items[0], int(key_items[1]), key_items[2]
        n_met, n_cov, rmet = freqinfo[fkey]
        wstr = "\t".join([chrom, str(pos), str(pos+1), ".", str(n_cov),
                          strand, str(pos), str(pos+1), "0,0,0", str(n_cov),
                          str(round(rmet*100))]) + "\n"
        # round(rmet*100) will cause imprecise, leads to lower performance,
        # case1: rmet=0.4652, case2: round(rmet*100)/100 turns to 0.47
        # e.g. rice1-1 CHH pearson cors with BS: case1 0.72196 vs case2 0.7219
        wf.write(wstr)
    wf.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bedfile", action="append", type=str, required=True,
                        help="methylbed file")
    parser.add_argument("--wfile", type=str, required=True,
                        help="")

    args = parser.parse_args()

    freqinfo = _get_combined_freq_info(args.bedfile)
    _write_freqinfo(freqinfo, args.wfile)


if __name__ == '__main__':
    main()
