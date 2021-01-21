#!/usr/bin/python
import argparse
import os

sep = "||"


# def read_mummer_repeat_file(mrfile):
#     region2group = {}
#     group2edges = {}
#     regions = set()
#     group_number = 0
#     with open(mrfile, "r") as rf:
#         for line in rf:
#             words = line.strip().split("\t")
#             if len(words) > 3:
#                 region1 = "\t".join(words[0:3])
#                 region2 = "\t".join(words[3:6])
#                 if region1 not in regions and region2 not in regions:
#                     regions.add(region1)
#                     regions.add(region2)
#
#                     gnumber = group_number
#                     region2group[region1] = gnumber
#                     region2group[region2] = gnumber
#                     group2edges[gnumber] = 1
#
#                     group_number += 1
#                 elif region1 in regions:
#                     regions.add(region2)
#
#                     gnumber = region2group[region1]
#                     region2group[region2] = gnumber
#                     group2edges[gnumber] += 1
#                 elif region2 in regions:
#                     regions.add(region1)
#
#                     gnumber = region2group[region2]
#                     region2group[region1] = gnumber
#                     group2edges[gnumber] += 1
#             else:
#                 region1 = "\t".join(words)
#                 if region1 not in regions:
#                     regions.add(region1)
#                 gnumber = group_number
#                 region2group[region1] = gnumber
#                 group2edges[gnumber] = 0
#
#                 group_number += 1
#     return region2group, group2edges
#
#
# def check_region_groups(region2group, group2edges):
#     group2regions = {}
#     for rgn in region2group.keys():
#         gnumber = region2group[rgn]
#         if gnumber not in group2regions.keys():
#             group2regions[gnumber] = set()
#         group2regions[gnumber].add(rgn)
#
#     for gnumber in group2regions.keys():
#         print(gnumber, len(group2regions[gnumber]), group2edges[gnumber])


def reformat_mummer_repeat_file(args):
    fname, fext = os.path.splitext(args.repeat_file)
    wfile = fname + ".reformat" + fext
    wf = open(wfile, "w")
    with open(args.repeat_file, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            if len(words) <= 3:
                continue
            region1 = sep.join(words[0:3])
            region2 = sep.join(words[3:6])
            wf.write("\t".join([region1, region2]) + "\n")
    wf.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat_file", type=str, required=True, help="from s1_extract_repeats_from_mummer.py, "
                                                                       "ref2ref file")
    args = parser.parse_args()

    # region2group, group2edges = read_mummer_repeat_file(args.repeat_file)
    # check_region_groups(region2group, group2edges)

    reformat_mummer_repeat_file(args)


if __name__ == '__main__':
    main()
