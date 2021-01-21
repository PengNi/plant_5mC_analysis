#! /usr/bin/python
import argparse

alphabeta = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def _read_one_nano_file(freqfile, cov_cf):
    chrom2poses = {}
    with open(freqfile, "r") as rf:
        # next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom = words[0]
            m_key = "\t".join([words[0], words[1], words[2]])
            cov = int(words[8])
            if cov >= cov_cf:
                if chrom not in chrom2poses.keys():
                    chrom2poses[chrom] = set()
                chrom2poses[chrom].add(m_key)
    return chrom2poses


def stat_poses(chrom2poses, contig_prefix, contig_names):
    poses_genome = set()
    totalcount = 0

    contigset = set(contig_names.strip().split(",")) if contig_names is not None else None

    for chrom in chrom2poses.keys():
        poses = chrom2poses[chrom]
        if contig_prefix is not None:
            if str(chrom).startswith(contig_prefix):
                poses_genome.update(poses)
        elif contig_names is not None:
            # contigset = contig_names
            if chrom in contigset:
                poses_genome.update(poses)
        else:
            poses_genome.update(poses)
        print("{}\t{}".format(len(poses), chrom))
        totalcount += len(poses)
    print("{}\tgenome(main)".format(len(poses_genome)))
    print("{}\tgenome(all)".format(totalcount))
    return poses_genome


def write_union_poses(poses, wfile):
    poses = sorted(list(poses))
    wf = open(wfile, "w")
    for mpos in poses:
        wf.write(mpos + "\n")
    wf.close()


def main():
    # python stats_count_c_of_nanopore.py --freqfile ~/tools/data/rice/shuidao1-1.guppy.pass.part2.CHH.bn13_sn16.arabnrice2-1.denoise_signal_bilstm.both_bilstm.50x_12345.freq.tsv --cov_cf 5 --contig_names 1,2,3,4,5,6,7,8,9,10,11,12 --wfile shuidao1-1.guppy.pass.part2.CHH.bn13_sn16.arabnrice2-1.denoise_signal_bilstm.both_bilstm.50x_12345.c_count.txt > shuidao1-1.guppy.pass.part2.CHH.bn13_sn16.arabnrice2-1.denoise_signal_bilstm.both_bilstm.50x_12345.c_count.log &
    parser = argparse.ArgumentParser("count motif num of each contig and genome in bisulfite sequencing")
    parser.add_argument("--freqfile", type=str, action="store", required=True, help="")
    parser.add_argument("--cov_cf", type=int, required=False, default=5, help="")
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12
    parser.add_argument("--wfile", type=str, required=True, default=None)

    args = parser.parse_args()
    print("=============coverage cutoff: 1")
    chrom2poses = _read_one_nano_file(args.freqfile, 1)
    poses = stat_poses(chrom2poses, args.contig_prefix, args.contig_names)
    print("=============coverage cutoff: {}".format(args.cov_cf))
    chrom2poses = _read_one_nano_file(args.freqfile, args.cov_cf)
    poses = stat_poses(chrom2poses, args.contig_prefix, args.contig_names)
    write_union_poses(poses, args.wfile)


if __name__ == '__main__':
    main()
