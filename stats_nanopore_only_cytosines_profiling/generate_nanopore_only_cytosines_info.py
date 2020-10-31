#! /usr/bin/python
import argparse
import os


def read_one_posfile(posfile):
    posinfo = dict()
    with open(posfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            posinfo["\t".join([words[0], words[1], words[2]])] = words[3]
    return posinfo


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


def add_info_to_posfile(args):
    freqfiles = args.nano_file
    freqinfo = dict()
    for freqfile in freqfiles:
        finfo_tmp = _read_one_mod_freq_file(freqfile)
        freqinfo.update(finfo_tmp)

    posinfo = read_one_posfile(args.posfile)
    for pos in list(posinfo.keys()):
        coverage, rmet = freqinfo[pos][5], freqinfo[pos][6]
        motif = posinfo[pos]
        posinfo[pos] = [coverage, rmet, motif]
    return posinfo


def write_posinfo(wfile, posinfo):
    wf = open(wfile, "w")
    pkeys = sorted(list(posinfo.keys()))
    for pkey in pkeys:
        wf.write("\t".join([str(pkey), ] + list(map(str, posinfo[pkey]))) + "\n")
    wf.close()


def main():
    parser = argparse.ArgumentParser("")
    parser.add_argument("--posfile", type=str, required=True,
                        help="from compare_detected_cytosines.py")
    parser.add_argument("--nano_file", type=str, action="append",
                        required=True, help="nanopore call_mods freq file")
    parser.add_argument("--wfile", type=str, required=True, default=None)

    args = parser.parse_args()
    posinfo = add_info_to_posfile(args)
    write_posinfo(args.wfile, posinfo)


if __name__ == '__main__':
    main()
