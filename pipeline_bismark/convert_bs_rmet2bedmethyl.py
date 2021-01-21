#! /usr/bin/python
import argparse


chromname_map_arab = {"NC_003070.9": "Chr1",
                      "NC_003071.7": "Chr2",
                      "NC_003074.8": "Chr3",
                      "NC_003075.7": "Chr4",
                      "NC_003076.8": "Chr5",
                      "NC_037304.1": "ChrM",
                      "NC_000932.1": "ChrPltd"}


def _read_one_bs_rmet_file(bs_file, args):
    freqinfo = {}
    with open(bs_file, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, coverage, rmet = words[0], int(words[1]), words[2], int(words[6]), float(words[7])
            if args.conv_chrom:
                chrom = chromname_map_arab[chrom]
            # mkey = "\t".join([chrom, str(pos), strand])
            mkey = (chrom, pos, strand)
            if coverage >= args.covcf:
                freqinfo[mkey] = (coverage, rmet)
    return freqinfo


def combine_bs_rmet_files(args):
    freqinfo_bs = []
    for bs_file in args.bsfile:
        freqtmp = _read_one_bs_rmet_file(bs_file, args)
        freqinfo_bs.append(freqtmp)
    freqinfo_bs_comb = dict()
    keyset = set()
    for freqinfotmp in freqinfo_bs:
        for fkey in freqinfotmp.keys():
            if fkey not in keyset:
                keyset.add(fkey)
                freqinfo_bs_comb[fkey] = []
            freqinfo_bs_comb[fkey].append(freqinfotmp[fkey])
    if args.sort:
        keyset = sorted(list(keyset))
    wf = open(args.wfile, "w")
    for fkey in keyset:
        infos = freqinfo_bs_comb[fkey]
        coverages, rmets = [], []
        for i in range(len(infos)):
            coverages.append(infos[i][0])
            rmets.append(infos[i][1])
        # chrom, pos, strand = fkey.split("\t")
        chrom, pos, strand = fkey
        pos = int(pos)
        # m_cov, m_rmet = int(round(sum(coverages)/float(len(coverages)), 0)), sum(rmets) / len(rmets)
        m_cov, m_rmet = sum(coverages), sum(rmets) / len(rmets)
        wf.write("\t".join([chrom, str(pos), str(pos+1), ".", str(m_cov), strand,
                            str(pos), str(pos + 1), "0,0,0", str(m_cov), str(int(round(m_rmet*100, 0)))]) + "\n")

    wf.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bsfile", type=str, action="append", required=True,
                        help="bs rmet file")
    parser.add_argument("--covcf", type=int, required=False, default=1,
                        help="")
    parser.add_argument("--wfile", type=str, required=True)
    parser.add_argument("--conv_chrom", action='store_true', default=False,
                        help="convert chrom names of arab")
    parser.add_argument('--sort', action='store_true', default=False, help="sort items in the result")

    args = parser.parse_args()
    combine_bs_rmet_files(args)


if __name__ == '__main__':
    main()
