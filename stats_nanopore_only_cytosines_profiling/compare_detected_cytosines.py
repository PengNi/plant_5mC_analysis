#! /usr/bin/python
import argparse
import os


def read_one_posfile(posfile):
    poses = set()
    with open(posfile, "r") as rf:
        for line in rf:
            poses.add(line.strip())
    return poses


def cmp_two_possets(posset1, posset2):
    posset1 = set(posset1)
    posset2 = set(posset2)
    num_inter = len(posset1.intersection(posset2))
    num_diff1 = len(posset1.difference(posset2))
    pos_diff2 = posset2.difference(posset1)
    num_diff2 = len(pos_diff2)
    return len(posset1), len(posset2), num_inter, num_diff1, num_diff2, pos_diff2


def cmp_detected_cytosines(args):
    wf = open(args.wfile, "w")
    print("==bs poses:")
    for bs_file in args.bs_file:
        print(bs_file)
    print("==nano poses:")
    for nanofile in args.nano_file:
        print(nanofile)
    print("\n==cmp stats==")
    motifs = str(args.motifs).split(",")
    combine_nums = [0, 0, 0, 0, 0]
    for m_idx in range(0, len(motifs)):
        bsfile = args.bs_file[m_idx]
        nanofile = args.nano_file[m_idx]
        print("{} vs {}:".format(bsfile, nanofile))
        len_1, len_2, len_inter, len_diff1, len_diff2, pos_diff2 = cmp_two_possets(read_one_posfile(bsfile),
                                                                                   read_one_posfile(nanofile))
        print("pos1: {}, pos2: {}, inter: {}, pos1.diff: {}, pos2.diff: {}".format(len_1, len_2,
                                                                                   len_inter,
                                                                                   len_diff1,
                                                                                   len_diff2))
        one_nums = [len_1, len_2, len_inter, len_diff1, len_diff2]
        combine_nums = [combine_nums[i] + one_nums[i] for i in range(0, len(combine_nums))]
        for d2pos in pos_diff2:
            wf.write(d2pos + "\t" + motifs[m_idx] + "\n")
    print("\n==combine: pos1: {}, pos2: {}, inter: {}, pos1.diff: {}, pos2.diff: {}".format(combine_nums[0],
                                                                                            combine_nums[1],
                                                                                            combine_nums[2],
                                                                                            combine_nums[3],
                                                                                            combine_nums[4],))
    wf.close()


def main():
    # python compare_detected_cytosines.py --bs_file ../stats_c_count/ninanjie-2.bs_3replicates.CG.main_genome.stat.txt --bs_file ../stats_c_count/ninanjie-2.bs_3replicates.CHG.main_genome.stat.txt --bs_file ../stats_c_count/ninanjie-2.bs_3replicates.CHH.main_genome.stat.txt --nano_file ../stats_c_count/ninanjie.nanopore50x.CG.main_genome.stat.txt --nano_file ../stats_c_count/ninanjie.nanopore50x.CHG.main_genome.stat.txt --nano_file ../stats_c_count/ninanjie.nanopore50x.CHH.main_genome.stat.txt --wfile ninanjie.c_count.bs_3reps_vs_nano.only_nanopore_detected_Cs.pos.txt > ninanjie.c_count.bs_3reps_vs_nano50x.log &
    # python compare_detected_cytosines.py --bs_file ../stats_c_count/ninanjie-2.bs_3replicates.CG.main_genome.stat.txt --bs_file ../stats_c_count/ninanjie-2.bs_3replicates.CHG.main_genome.stat.txt --bs_file ../stats_c_count/ninanjie-2.bs_3replicates.CHH.main_genome.stat.txt --nano_file ../stats_c_count/athaliana.guppy.pass.part2.CG.dp2_p0.8_50x12345.main_genome.stat.txt --nano_file ../stats_c_count/athaliana.guppy.pass.part2.CHG.dp2_p0.8_50x12345.main_genome.stat.txt --nano_file ../stats_c_count/athaliana.guppy.pass.part2.CHH.dp2_p0.8_50x12345.main_genome.stat.txt --wfile ninanjie.c_count.bs_3reps_vs_nano50x_dp2_p0.8.only_nanopore_detected_Cs.pos.txt > ninanjie.c_count.bs_3reps_vs_nano50x_dp2_p0.8.log &
    parser = argparse.ArgumentParser("must add bs_file/nano_file in the same order, from stats_c_count")
    parser.add_argument("--bs_file", type=str, action="append",
                        required=True, help="bs detected poses")
    parser.add_argument("--nano_file", type=str, action="append",
                        required=True, help="nanopore detected poses")
    parser.add_argument("--motifs", type=str, default="CG,CHG,CHH",
                        required=False, help="")
    parser.add_argument("--wfile", type=str, required=True, default=None)

    args = parser.parse_args()
    cmp_detected_cytosines(args)


if __name__ == '__main__':
    main()
