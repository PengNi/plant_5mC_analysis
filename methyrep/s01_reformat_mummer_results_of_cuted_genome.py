#!/usr/bin/python
import argparse
import os


class MumRepeat:
    def __init__(self, fields):
        self._refstart = int(fields[0])
        self._refend = int(fields[1])
        self._querystart = int(fields[2])
        self._queryend = int(fields[3])
        self._ref_regionlen = int(fields[4])
        self._query_regionlen = int(fields[5])
        self._identity = float(fields[6])
        self._refid = fields[11]
        self._queryid = fields[12]


def change_chr_coords(words):
    ref_s, ref_e = int(words[0]), int(words[1])
    query_s, query_e = int(words[2]), int(words[3])
    refid = words[11]
    queryid = words[12]

    ref_wrds = refid.split("_")
    refname, refstart, refend = "_".join(ref_wrds[:-2]), int(ref_wrds[-2]), int(ref_wrds[-1])
    query_wrds = queryid.split("_")
    queryname, querystart, queryend = "_".join(query_wrds[:-2]), int(query_wrds[-2]), int(query_wrds[-1])

    ref_s += refstart
    ref_e += refstart
    query_s += querystart
    query_e += querystart
    return [str(ref_s), str(ref_e), str(query_s), str(query_e)] + words[4:11] + [refname, queryname]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--crds_ref_file", type=str, required=True, help="nucmer show-coords ref2ref file, "
                                                                         "or contig2contig file")

    args = parser.parse_args()
    fname, fext = os.path.splitext(args.crds_ref_file)
    wfile = fname + ".chr_offset2absloc" + fext
    wf = open(wfile, "w")
    with open(args.crds_ref_file, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            if len(words) < 13:
                wf.write(line)
            else:
                words_r = change_chr_coords(words)
                wf.write("\t".join(words_r) + "\n")
    wf.close()


if __name__ == '__main__':
    main()
