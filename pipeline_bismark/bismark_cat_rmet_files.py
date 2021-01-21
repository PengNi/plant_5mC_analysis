#! /usr/bin/python
import argparse


def _read_one_bs_rmet_file(bs_file):
    freqinfo = {}
    with open(bs_file, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, coverage, met = words[0], words[1], words[2], int(words[6]), int(words[4])
            mkey = "\t".join([chrom, str(pos), strand])
            # if coverage >= args.covcf:
            freqinfo[mkey] = (coverage, met)
    return freqinfo


def combine_bs_rmet_files(args):
    freqinfo_bs = []
    for bs_file in args.bsfile:
        freqtmp = _read_one_bs_rmet_file(bs_file)
        freqinfo_bs.append(freqtmp)
    freqinfo_bs_comb = dict()
    keyset = set()
    for freqinfotmp in freqinfo_bs:
        for fkey in freqinfotmp.keys():
            if fkey not in keyset:
                keyset.add(fkey)
                freqinfo_bs_comb[fkey] = []
            freqinfo_bs_comb[fkey].append(freqinfotmp[fkey])
    wf = open(args.wbed, "w")
    c_cnts = 0
    for fkey in freqinfo_bs_comb.keys():
        infos = freqinfo_bs_comb[fkey]
        coverages, mets = [], []
        for i in range(len(infos)):
            coverages.append(infos[i][0])
            mets.append(infos[i][1])
        chrom, pos, strand = fkey.split("\t")
        pos = int(pos)
        # m_cov, m_met = sum(coverages), sum(mets)
        # m_rmet = float(m_met)/m_cov
        m_cov = 1
        for cov in coverages:
            if cov >= args.covcf:
                m_cov = args.covcf
                c_cnts += 1
                break
        m_rmet = 0
        wf.write("\t".join([chrom, str(pos), str(pos+1), ".", str(m_cov), strand,
                            str(pos), str(pos + 1), "0,0,0", str(m_cov), str(int(round(m_rmet*100, 0)))]) + "\n")

    wf.close()
    print("{} sites with at least {} reads".format(c_cnts, args.covcf))


def main():
    parser = argparse.ArgumentParser("concat coverage with coverage_cutoff, "
                                     "count sites with  at least covcf reads in at least 1 bsfile")
    parser.add_argument("--bsfile", type=str, action="append", required=True,
                        help="bs rmet file")
    parser.add_argument("--covcf", type=int, required=False, default=5,
                        help="")
    parser.add_argument("--wbed", type=str, required=True)

    args = parser.parse_args()
    combine_bs_rmet_files(args)


if __name__ == '__main__':
    main()
