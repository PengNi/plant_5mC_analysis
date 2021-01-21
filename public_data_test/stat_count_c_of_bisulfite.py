#! /usr/bin/python
import argparse

alphabeta = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


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
def read_bsmap_file(bsmapfile):
    methylinfo = dict()
    with open(bsmapfile, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, motif, rmet, effcov = words[0], int(words[1]) - 1, words[2], words[3], \
                float(words[4]), float(words[5])
            # m_cov, all_cov = int(words[6]), int(words[7])
            # methylinfo[motif].append([chrom, pos, strand, effcov, rmet, m_cov, all_cov])
            if chrom not in methylinfo.keys():
                methylinfo[chrom] = {}
            if motif not in methylinfo[chrom].keys():
                methylinfo[chrom][motif] = []
            methylinfo[chrom][motif].append((pos, strand, effcov))
    return methylinfo


def stat_poses(methylinfo, contig_names, cov_cf):
    contigset = set(contig_names.strip().split(",")) if contig_names is not None else None

    print("=====cov_cf: {}".format(cov_cf))
    coverages = []
    total_motifs = [0, 0, 0]
    main_motifs = [0, 0, 0]
    motifs = ["CG", "CHG", "CHH"]
    print("chrom\t" + "\t".join(motifs))
    for chrom in methylinfo.keys():
        chrom_cnt = [0, 0, 0]
        for mfidx in range(len(motifs)):
            if motifs[mfidx] in methylinfo[chrom].keys():
                chrom_motif_cnt = 0
                for mitem in methylinfo[chrom][motifs[mfidx]]:
                    pos, strand, effcov = mitem
                    coverages.append(effcov)
                    if effcov >= cov_cf:
                        chrom_motif_cnt += 1
                chrom_cnt[mfidx] = chrom_motif_cnt

                total_motifs[mfidx] += chrom_motif_cnt
                if contigset is not None and chrom in contigset:
                    main_motifs[mfidx] += chrom_motif_cnt
        print("{}\t".format(chrom) + "\t".join(list(map(str, chrom_cnt))))
    print("genome(main)\t" + "\t".join(list(map(str, main_motifs))))
    print("genome(all)\t" + "\t".join(list(map(str, total_motifs))))
    print("====all Cs mean coverage: {}\n".format(float(sum(coverages))/len(coverages)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bsmapfile", type=str, action="store", required=True,
                        help="bsmapfile")
    parser.add_argument("--cov_cf", type=int, required=False, default=5, help="")
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12

    args = parser.parse_args()
    methylinfo = read_bsmap_file(args.bsmapfile)
    stat_poses(methylinfo, args.contig_names, args.cov_cf)
    stat_poses(methylinfo, args.contig_names, 1)


if __name__ == '__main__':
    main()
