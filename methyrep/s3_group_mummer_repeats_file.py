#!/usr/bin/python
import argparse
import networkx as nx
import os

sep = "||"


def str2bool(v):
    # susendberg's function
    return v.lower() in ("yes", "true", "t", "1")


def remove_replicate_nodes(nodes, olap_fac=0.0, is_combine=False):
    nodes = sorted(list(nodes))
    uni_nodes = [nodes[0], ]
    for i in range(1, len(nodes)):
        last_chrom, last_s, last_e = uni_nodes[-1]

        curr_chrom, curr_s, curr_e = nodes[i]

        if curr_chrom == last_chrom:
            dist_cf = int(0.1 * abs(last_e - last_s))
            if abs(curr_s - last_s) <= dist_cf and abs(curr_e - last_e) <= dist_cf:
                if is_combine:
                    end_max = max(last_e, curr_e)
                    start_min = min(last_s, curr_s)
                    uni_nodes[-1] = (curr_chrom, start_min, end_max)
                else:
                    if abs(curr_e - curr_s) > abs(last_e - last_s):
                        uni_nodes[-1] = (curr_chrom, curr_s, curr_e)
            else:
                # remove overlap repeats
                end_min = min(last_e, curr_e)
                start_max = max(last_s, curr_s)
                if end_min - start_max < olap_fac * (last_e - last_s):
                    uni_nodes.append(nodes[i])
        else:
            uni_nodes.append(nodes[i])
    return uni_nodes


def group_repeats(repeat_file, is_combine_near=False):
    rgraph = nx.read_edgelist(repeat_file, delimiter="\t")

    fname, fext = os.path.splitext(repeat_file)
    group_dir = fname + ".group"
    if not os.path.exists(group_dir):
        os.mkdir(group_dir)
    else:
        for file in os.listdir(group_dir):
            os.remove(group_dir + "/" + file)
    rccidx = 1
    for rcc in sorted(nx.connected_components(rgraph), key=len, reverse=True):
        subrg = rgraph.subgraph(rcc).copy()

        words = str(list(rcc)[0]).split(sep)
        _, start, end = words[0], int(words[1]), int(words[2])
        region_len = end - start + 1

        srgnodes = []
        for srgn in rcc:
            words = str(srgn).split(sep)
            chrom, start, end = words[0], int(words[1]), int(words[2])
            srgnodes.append((chrom, start, end))
        srgnodes = sorted(list(srgnodes))
        srgnodes = remove_replicate_nodes(srgnodes, is_combine_near)
        if len(srgnodes) > 1:
            wf = open(group_dir + "/" + str(rccidx).zfill(4) + "_n" + str(len(srgnodes)) +
                      "_region" + str(region_len) + ".nodes.tsv", "w")
            for srge in subrg.edges:
                wf.write("##" + str(srge) + "\n")
            for srgn in sorted(srgnodes):
                wf.write("\t".join(list(map(str, srgn))) + "\n")
            wf.close()
            rccidx += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat_file", type=str, required=True, help="from s2_reformat_mummer_repeats.py")
    parser.add_argument("--is_comb_near", type=str, required=False, default="no",
                        help="is combine overlap repeats. if no, only keep the 1st one")
    parser.add_argument("--olap_fac", type=float, required=False, default=0.0,
                        help="0-1, 0.0 means no overlap for two repeats allowed")

    args = parser.parse_args()

    group_repeats(args.repeat_file, str2bool(args.is_comb_near))


if __name__ == '__main__':
    main()
