#! /usr/bin/python
import argparse
import os

key_sep = "||"


def _bin(key, key2bsbin):
    return key2bsbin[key]


def _read_keys(bsbinfile):
    key2bin = dict()
    with open(bsbinfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            key2bin[words[0]] = words[2]
    return key2bin


def _read_nano_file(nanofile, cov_cf, bskey2bin):
    bskeys = set(bskey2bin.keys())
    items = []
    with open(nanofile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            if str(nanofile).endswith(".bed"):
                chrom, pos, coverage, rmet = words[0], words[1], float(words[9]), float(words[10]) / 100
            else:
                chrom, pos, coverage, rmet = words[0], words[1], int(words[8]), float(words[9])
            if coverage < cov_cf:
                continue
            keynano = key_sep.join([chrom, pos])
            if keynano not in bskeys:
                continue
            site_bin = _bin(keynano, bskey2bin)
            items.append((keynano, rmet, site_bin))
    return items


def write_items(items, wfile, technique="bisulfite"):
    wf = open(wfile, "w")
    for it in items:
        wf.write("\t".join(list(map(str, it)) + [technique, ]) + "\n")
    wf.flush()
    wf.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bsbin", type=str, action="store", required=True,
                        help="bsbin file , from stat_methyfreqs_in_bins.py")
    parser.add_argument("--nanofile", type=str, action="store", required=True,
                        help="")
    parser.add_argument("--cov_cf", type=int, default=5, required=False,
                        help="")
    parser.add_argument("--techname", type=str, action="store", required=True,
                        help="")
    parser.add_argument("--wfile", type=str, required=False, help="")

    args = parser.parse_args()

    bskeys = _read_keys(args.bsbin)
    nanoitems = _read_nano_file(args.nanofile, args.cov_cf, bskeys)
    print("file: {}, covcf: {}, cmp_bsfile: {}\n\titems: {}".format(args.nanofile,
                                                                    args.cov_cf,
                                                                    args.bsbin,
                                                                    len(nanoitems)))
    write_items(nanoitems, args.wfile, args.techname)


if __name__ == '__main__':
    main()
