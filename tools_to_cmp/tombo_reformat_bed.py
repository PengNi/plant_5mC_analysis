#! /usr/bin/python
import argparse


def read_coverage_bedgraph(covbed):
    pos2coverage = dict()
    with open(covbed, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom, start, end, coverage = words[0], int(words[1]), int(words[2]), int(words[3])
            for pos in range(start, end):
                poskey = "\t".join([chrom, str(pos)])
                pos2coverage[poskey] = coverage
    return pos2coverage


def read_rmet_bed(rmetbed):
    pos2rmet = dict()
    with open(rmetbed, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            chrom, start, end, rmet = words[0], int(words[1]), int(words[2]), float(words[4])
            poskey = "\t".join([chrom, str(start)])
            pos2rmet[poskey] = rmet
    return pos2rmet


def combine_cov_rmet(pos2rmet, pos2cov, strand):
    bedeles = []
    for poskey in pos2rmet.keys():
        chrom = poskey.split("\t")[0]
        start = int(poskey.split("\t")[1])
        rmet = pos2rmet[poskey]
        rmet_percent = int(round(rmet*100, 0))
        coverage = pos2cov[poskey]
        bedeles.append((chrom, start, start+1, ".", coverage, strand,
                        start, start+1, "0,0,0", coverage, rmet_percent))
    return bedeles


def write_bedeles_to_file(bedeles, wfile):
    bedeles = sorted(bedeles)
    wf = open(wfile, "w")
    for bedele in bedeles:
        wf.write("\t".join(list(map(str, list(bedele)))) + "\n")
    wf.close()


def main():
    parser = argparse.ArgumentParser("combine tombo.cov.bed and tombo.rmet.bed to bedmethyl")
    parser.add_argument("--cov_plus", type=str, required=True, help="coverage bedgraph plus")
    parser.add_argument("--cov_minus", type=str, required=True, help="coverage bedgraph minus")
    parser.add_argument("--rmet_plus", type=str, required=True, help="rmet bed plus")
    parser.add_argument("--rmet_minus", type=str, required=True, help="rmet bed minus")
    parser.add_argument("--wfile", type=str, required=True, help="write file path, bedmethyl")

    args = parser.parse_args()
    print("raed plus strand ===")
    bedeles_plus = combine_cov_rmet(read_rmet_bed(args.rmet_plus),
                                    read_coverage_bedgraph(args.cov_plus), "+")
    print("raed minus strand ===")
    bedeles_minus = combine_cov_rmet(read_rmet_bed(args.rmet_minus),
                                     read_coverage_bedgraph(args.cov_minus), "-")
    print("combine to bed ===")
    bedeles = bedeles_plus + bedeles_minus
    write_bedeles_to_file(bedeles, args.wfile)


if __name__ == '__main__':
    main()
