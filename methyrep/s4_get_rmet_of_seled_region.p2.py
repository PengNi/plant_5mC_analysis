#!/usr/bin/python
import argparse
import os


chromname_map_arab = {"NC_003070.9": "Chr1",
                      "NC_003071.7": "Chr2",
                      "NC_003074.8": "Chr3",
                      "NC_003075.7": "Chr4",
                      "NC_003076.8": "Chr5",
                      "NC_037304.1": "ChrM",
                      "NC_000932.1": "ChrC"}
chromname_map_arab2 = {"Chr1": "NC_003070.9",
                       "Chr2": "NC_003071.7",
                       "Chr3": "NC_003074.8",
                       "Chr4": "NC_003075.7",
                       "Chr5": "NC_003076.8",
                       "ChrM": "NC_037304.1",
                       "ChrC": "NC_000932.1"}


def str2bool(v):
    # susendberg's function
    return v.lower() in ("yes", "true", "t", "1")


def read_regions(region_file):
    region2abslocs = dict()
    with open(region_file, "r") as rf:
        for line in rf:
            if line.startswith("##"):
                continue
            words = line.strip().split("\t")
            region2abslocs[words[0]] = list(map(int, words[1:]))
    return region2abslocs


def read_bs_rmet_info(bs_rmet_file):
    # chromosome	pos	strand	pos_in_strand	met	unmet	coverage	Rmet	cRmet	p.value	q.value	Rtype
    posinfo = []
    with open(bs_rmet_file, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            try:
                chrom, pos, strand, rmet, coverage = words[0], int(words[1]), words[2], float(words[7]), int(words[6])
                posinfo.append((chrom, pos, strand, rmet, coverage))
            except ValueError:
                continue
    return posinfo


def read_bs_bed_file(bs_bed):
    # chrom, pos, pos_end, _, coverage, strand, _, _, _, coverage, rmet*100
    posinfo = []
    with open(bs_bed, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, rmet, coverage = words[0], int(words[1]), words[5], \
                float(words[10]) / 100, int(words[9])
            if chrom in chromname_map_arab2.keys():
                chrom = chromname_map_arab2[chrom]
            posinfo.append((chrom, pos, strand, rmet, coverage))
    return posinfo


def read_deepsignal_rmet_file(nano_rmet_file):
    posinfo = []
    with open(nano_rmet_file, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, rmet, coverage = words[0], int(words[1]), words[2], \
                float(words[9]), int(words[8])
            posinfo.append((chrom, pos, strand, rmet, coverage))
    return posinfo


def read_rmet_info(rmet_file, rmet_type):
    if rmet_type == "bs_rmet":
        return read_bs_rmet_info(rmet_file)
    elif rmet_type == "bedmethyl":
        return read_bs_bed_file(rmet_file)
    elif rmet_type == "deepsignal":
        return read_deepsignal_rmet_file(rmet_file)
    else:
        return None


def posinfo_list2dict(posinfo):
    pos2rmetinfo = {}
    poses = set()
    for pinfo in posinfo:
        chrom, pos, strand, rmet, coverage = pinfo
        pos2rmetinfo[(chrom, pos)] = (rmet, coverage)
        poses.add((chrom, pos))
    return pos2rmetinfo, poses


def get_region_poses(regions, pos2rmetinfo, poses):
    regionkeys = sorted(list(regions.keys()))
    region2rmetstr = ""
    # head str ====================================================
    region_eles = regionkeys[-1].split("_")
    rchrom, rstart, rend, rstrand = "_".join(region_eles[:-3]), int(region_eles[-3]), \
                                    int(region_eles[-2]), region_eles[-1]
    rlocs = regions[regionkeys[-1]]
    offset_locs = []
    if len(rlocs) == 1:
        offset_locs.append(rlocs[0] - (rstart - 1))
    else:
        if rlocs[0] < rlocs[1]:
            offset_locs = [x - (rstart - 1) for x in rlocs]
        else:
            offset_locs = [(rend - 1) - x for x in rlocs]
    region2rmetstr += "region\t{}\n".format("\t".join(list(map(str, offset_locs))))
    for region in regionkeys:
        region_eles = region.split("_")
        rchrom = "_".join(region_eles[:-3])
        region2rmetstr += region + "\t"
        rrmets = []
        for rloc in regions[region]:
            if (rchrom, rloc) in poses:
                rmettmp, covertmp = pos2rmetinfo[(rchrom, rloc)]
                # if covertmp >= 5:
                #     rrmets.append(rmettmp)
                # else:
                #     rrmets.append(-0.01)
                rrmets.append(",".join([str(rmettmp), str(covertmp)]))
            else:
                rrmets.append(",".join([str(-0.01), str(0)]))
        region2rmetstr += "\t".join(list(map(str, rrmets))) + "\n"
    return region2rmetstr


def get_bs_rmet_of_seled_region(args):
    posinfo = read_rmet_info(args.rmet, args.rmet_type)
    pos2rmetinfo, poses = posinfo_list2dict(posinfo)

    outdir = args.mum_region.rstrip("/") + "." + args.rmet_type + args.suffix
    # outdir = args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        for file in os.listdir(outdir):
            os.remove(outdir + "/" + file)
    for rfile in os.listdir(args.mum_region):
        fname, fext = os.path.splitext(rfile)
        regions = read_regions(args.mum_region + "/" + rfile)
        region2posinfo = get_region_poses(regions, pos2rmetinfo, poses)
        wfile = outdir + "/" + fname + ".repeat_region.txt"
        with open(wfile, "w") as wf:
            wf.write(region2posinfo)


def main():
    parser = argparse.ArgumentParser("input: bs rmet reformat file, mummer region file (chrom, start(1-based), "
                                     "end)")
    parser.add_argument("--rmet", type=str, required=True, help="")
    parser.add_argument("--rmet_type", type=str, required=False,
                        default="bs_rmet", choices=["bs_rmet", "bedmethyl", "deepsignal"], help="")
    parser.add_argument("--mum_region", type=str, required=True,
                        help="dir, from s4_get_rmet_of_seled_region.p1.py")
    parser.add_argument("--suffix", type=str, default="", required=False,
                        help="outdir suffix, rmet data type/sample info")
    # parser.add_argument("--outdir", type=str, required=True, help="")

    args = parser.parse_args()
    get_bs_rmet_of_seled_region(args)


if __name__ == '__main__':
    main()
