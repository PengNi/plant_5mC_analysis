#! /usr/bin/python
import argparse
import os

key_sep = "||"


def _bin(rmet):
    if 0 <= rmet <= 0.3:
        return "low"
    elif 0.3 < rmet < 0.7:
        return "mediate"
    elif 0.7 <= rmet <= 1:
        return "high"
    return "NA"


def _read_bs_bed(bsfile, cov_cf, contig_names):
    contigset = set(contig_names.strip().split(",")) if contig_names is not None else None

    bsitems = []
    with open(bsfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, coverage, rmet = words[0], words[1], float(words[9]), float(words[10])/100
            if contigset is not None and chrom not in contigset:
                continue
            if coverage < cov_cf:
                continue
            site_bin = _bin(rmet)
            bsitems.append((key_sep.join([chrom, pos]), rmet, site_bin))
    return bsitems


def write_items(items, wfile, technique="bisulfite"):
    wf = open(wfile, "w")
    for it in items:
        wf.write("\t".join(list(map(str, it)) + [technique, ]) + "\n")
    wf.flush()
    wf.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bsbed", type=str, action="store", required=True,
                        help="bsbedfile")
    parser.add_argument("--cov_cf_bs", type=int, default=10, required=False,
                        help="")
    # B1,B2,B3,B4,B5,B6,B7,B8  # 1,2,3,4,5,6,7,8,9,10,11,12
    parser.add_argument("--contig_names", type=str, required=False, default=None)
    parser.add_argument("--wfile", type=str, required=False, help="")

    args = parser.parse_args()

    bsitems = _read_bs_bed(args.bsbed, args.cov_cf_bs, args.contig_names)
    print("file: {}, covcf: {}, contig_names: {}\n\titems: {}".format(args.bsbed,
                                                                      args.cov_cf_bs,
                                                                      args.contig_names,
                                                                      len(bsitems)))
    write_items(bsitems, args.wfile)


if __name__ == '__main__':
    main()
