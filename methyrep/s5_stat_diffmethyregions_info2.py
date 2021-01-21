import argparse
import os

rmet_abs_cf = 0.5


def _parse_regionid(regionid):
    words = regionid.split("_")
    start, end, strand = int(words[-3]), int(words[-2]), words[-1]
    contig = "_".join(words[:-3])
    return contig, start, end, strand


def read_one_repeat_sites_info(repeat_file):
    posinfo = dict()
    posinfo["repeat_group"] = os.path.basename(repeat_file).split(".repeat_region.txt")[0]
    with open(repeat_file, "r") as rf:
        header = next(rf)
        words = header.strip().split("\t")
        offset_poses = [int(x) for x in words[1:]]
        rmets_list, regionids = [], []
        for line in rf:
            words = line.strip().split("\t")
            region_id = words[0]
            rmetncovs = [x.split(",") for x in words[1:]]
            rmets = []
            for rmetcov in rmetncovs:
                rmet, cov = float(rmetcov[0]), int(rmetcov[1])
                if cov >= 5:
                    rmets.append(rmet)
                else:
                    rmets.append(-0.01)
            rmets_list.append(rmets)
            regionids.append(region_id)
        posinfo["offset_poses"] = offset_poses
        posinfo["region_ids"] = regionids
        posinfo["rmets_list"] = rmets_list
    return posinfo


def read_region_locs(region_pos_fp):
    regionid2locs = dict()
    with open(region_pos_fp, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            regionid2locs[words[0]] = list(map(int, words[1:]))
    return regionid2locs


def parse_repeatgroupname(rgname):
    words = rgname.split(".")[0].split("_")
    rgnumber, groupsize, rglen = words[0], words[1], words[2]
    groupsize = int(groupsize[1:])
    rglen = int(rglen[6:])
    return groupsize, rglen


def stat_repeat_regions(posinfo, region_locs_info, abs_cf=0.5):
    repeat_group_id = posinfo["repeat_group"]
    offset_poses = posinfo["offset_poses"]
    region_ids = posinfo["region_ids"]
    rmets_list = posinfo["rmets_list"]

    diff_poses = []
    all2_poses = []
    cov_poses = []
    methylevel = None
    if len(region_ids) == 2:
        rmet1, rmet2 = rmets_list[0], rmets_list[1]
        methylevel = dict()
        methylevel1, methylevel2 = [-1] * len(rmet1), [-1] * len(rmet1)
        for i in range(len(rmet1)):
            all2_poses.append(i)
            if rmet1[i] > -0.01 and rmet2[i] > -0.01:
                cov_poses.append(i)
            if rmet1[i] > -0.01 and rmet2[i] > -0.01 and abs(rmet1[i] - rmet2[i]) >= abs_cf:
                diff_poses.append(i)
                if rmet1[i] > rmet2[i]:
                    methylevel1[i], methylevel2[i] = 1, 0
                elif rmet1[i] < rmet2[i]:
                    methylevel1[i], methylevel2[i] = 0, 1
        methylevel[0] = methylevel1
        methylevel[1] = methylevel2
    is_2regions = 1 if len(all2_poses) > 0 else 0
    is_2regions_1moreCs = 1 if len(all2_poses) >= 1 else 0
    cov_ratio = 0
    if is_2regions_1moreCs == 1:
        cov_ratio = len(cov_poses) / float(len(all2_poses))
    is_2regions_cmp = 1 if cov_ratio >= 0.5 else 0
    is_2regions_diff = 1 if is_2regions_cmp and len(diff_poses) / float(len(cov_poses)) >= 0.1 else 0
    is_2regions_diff_sig = 1 if is_2regions_cmp and len(diff_poses) / float(len(cov_poses)) >= 0.5 else 0

    poses_all = []
    poses_all2 = []
    poses_diff = []
    # chrom, pos, strand, repeatid, repeatgroupid
    for i in range(0, len(region_ids)):
        rgid = region_ids[i]
        regionlocs = region_locs_info[rgid]
        chrom, start, end, strand = _parse_regionid(rgid)

        pos_abs = [regionlocs[x] for x in range(len(offset_poses))]
        poses_all += [(chrom, p, strand, rgid, repeat_group_id, -1) for p in pos_abs]

        if len(all2_poses) > 0:
            # pos_abs2 = [regionlocs[x] for x in all2_poses]
            # methylevel_i = [methylevel[i][x] for x in all2_poses]
            poses_all2 += [(chrom, regionlocs[x], strand, rgid, repeat_group_id, methylevel[i][x]) for x in all2_poses]
        if len(diff_poses) > 0:
            # posdiff_abs = [regionlocs[x] for x in diff_poses]
            poses_diff += [(chrom, regionlocs[x], strand, rgid, repeat_group_id, methylevel[i][x]) for x in diff_poses]
    region_info = ""
    if is_2regions:
        groupsize, rglen = parse_repeatgroupname(repeat_group_id)
        if len(all2_poses) == 0:
            ratio_to_all = 0
            ratio_c2a = 0
        else:
            ratio_to_all = len(diff_poses)/float(len(all2_poses))
            ratio_c2a = len(cov_poses)/float(len(all2_poses))
        if len(cov_poses) == 0:
            ratio_to_cov = 0
        else:
            ratio_to_cov = len(diff_poses) / float(len(cov_poses))
        region_info = "\t".join([repeat_group_id, str(len(all2_poses)), str(len(cov_poses)),
                                 str(len(diff_poses)), str(ratio_c2a), str(ratio_to_all),
                                 str(ratio_to_cov),
                                 str(rglen)]) + "\n"
    return poses_all, poses_all2, poses_diff, is_2regions, is_2regions_1moreCs, \
        is_2regions_cmp, is_2regions_diff, is_2regions_diff_sig, region_info


# write results ==
def write_regions_to_file(regions, wfile):
    regions = sorted(list(regions))
    wf = open(wfile, "w")
    for rginfo in regions:
        wf.write("\t".join(list(rginfo)) + "\n")
    wf.close()


def write_regions_to_bed(regions, wfile, rg2methylevel):
    regions = sorted(list(regions))
    region_bed = []
    for rginfo in regions:
        groupid, regiontmp = rginfo
        chrom, start, end, strand = _parse_regionid(regiontmp)
        mlevel = rg2methylevel[(regiontmp, groupid)]
        region_bed.append((chrom, start-1, end, groupid, str(mlevel), strand))
    wf = open(wfile, "w")
    region_bed = sorted(region_bed)
    for rginfo in region_bed:
        wf.write("\t".join(list(map(str, rginfo))) + "\n")
    wf.close()


def write_poses_to_file(poses, wfile):
    wf = open(wfile, "w")
    for posinfo in poses:
        wf.write("\t".join(list(map(str, posinfo))) + "\n")
    wf.close()


def write_poses_to_file_bed(poses, wfile):
    wf = open(wfile, "w")
    poses = sorted(poses, key=lambda x: (x[0], x[1]))
    for posinfo in poses:
        chrom, p, strand, rgid, repeat_group_id, methylevel = posinfo
        wf.write("\t".join([chrom, str(p), str(p+1), rgid, str(methylevel), strand]) + "\n")
    wf.close()


def _comb_methylevels(methylevels):
    numlevel = dict()
    numlevel[0] = 0
    numlevel[1] = 0
    numlevel[-1] = 0
    for x in methylevels:
        numlevel[x] += 1
    if numlevel[-1] == len(methylevels):
        return -1
    return 1 if numlevel[1] > numlevel[0] else 0


def write_posregions_to_file_bed(poses, wfile):
    poses = sorted(poses, key=lambda x: (x[0], x[1]))
    regionnames = set()
    regionname2methylevel2 = dict()
    for posinfo in poses:
        chrom, p, strand, rgid, repeat_group_id, methylevel = posinfo
        if (rgid, repeat_group_id) not in regionnames:
            regionnames.add((rgid, repeat_group_id))
            regionname2methylevel2[(rgid, repeat_group_id)] = []
        regionname2methylevel2[(rgid, repeat_group_id)].append(methylevel)
    region_bed = []
    rg2methylevel = dict()
    for (rgid, repeat_group_id) in regionnames:
        contig, start, end, strand = _parse_regionid(rgid)
        methylevels = regionname2methylevel2[(rgid, repeat_group_id)]
        comb_mlevel = _comb_methylevels(methylevels)
        rg2methylevel[(rgid, repeat_group_id)] = comb_mlevel
        region_bed.append((contig, start-1, end, repeat_group_id, str(comb_mlevel), strand))
    region_bed = sorted(region_bed, key=lambda x: (x[0], x[1]))
    wf = open(wfile, "w")
    for rginfo in region_bed:
        wf.write("\t".join(list(map(str, rginfo))) + "\n")
    wf.close()
    return rg2methylevel


def filter_poses_by_regionid(poses, regionsinfo):
    regionids = set()
    for rinfo in regionsinfo:
        regionids.add(rinfo[1])
    poses_filt = []
    for posinfo in poses:
        if posinfo[3] in regionids:
            poses_filt.append(posinfo)
    return poses_filt


def main():
    parser = argparse.ArgumentParser("find diff methy sites from mum_repeats")
    parser.add_argument("--region_dir", type=str, required=True,
                        help="out dir from s4_get_rmet_of_seled_region.p2.py")
    parser.add_argument("--region_pos", type=str, required=True,
                        help="out dir from s4_get_rmet_of_seled_region.p1.py")
    parser.add_argument("--outdir", type=str, required=False, default=None, help="")
    parser.add_argument("--rmet_cf", type=float, required=False, default=0.5,
                        help="rmet abs cutoff to see if a site is diff_methyl in two repeats")

    args = parser.parse_args()

    if args.outdir is not None and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    cnt_100, cnt_1k, cnt_10k = 0, 0, 0
    all2_repeats_stats = ""
    allrepeats_poses_all, allrepeats_poses_all2, allrepeats_poses_diff = [], [], []
    cnt_all = 0
    cnt_diff = 0
    cnt_all2, cnt_all2_1mCs, cnt_all2_cmp, cnt_all2_diff, cnt_all2_diff_sig = 0, 0, 0, 0, 0
    regions_cmp, regions_diff, regions_diff_sig = [], [], []
    for file in os.listdir(args.region_dir):
        fname_prefix = os.path.basename(file).split(".repeat_region.txt")[0]

        groupsize, rglen = parse_repeatgroupname(fname_prefix)
        if groupsize == 2:
            if rglen >= 100:
                cnt_100 += 1
            if rglen >= 1000:
                cnt_1k += 1
            if rglen >= 10000:
                cnt_10k += 1

        regionloc_fp = args.region_pos + "/" + fname_prefix + ".txt"
        regionid2locs = read_region_locs(regionloc_fp)
        fpath = "/".join([args.region_dir, file])
        posinfo = read_one_repeat_sites_info(fpath)

        poses_all, poses_all2, poses_diff, is_2regions, is_2regions_1moreCs, \
           is_2regions_cmp, is_2regions_diff, is_2regions_diff_sig, region_info = stat_repeat_regions(posinfo,
                                                                                                      regionid2locs,
                                                                                                      args.rmet_cf)
        all2_repeats_stats += region_info

        allrepeats_poses_all += poses_all
        allrepeats_poses_all2 += poses_all2
        allrepeats_poses_diff += poses_diff

        cnt_diff = cnt_diff + 1 if len(poses_diff) > 0 else cnt_diff

        cnt_all += 1
        cnt_all2 += is_2regions
        cnt_all2_1mCs += is_2regions_1moreCs
        cnt_all2_cmp += is_2regions_cmp
        cnt_all2_diff += is_2regions_diff
        cnt_all2_diff_sig += is_2regions_diff_sig

        if is_2regions_cmp:
            for regionid in posinfo["region_ids"]:
                regions_cmp.append((fname_prefix, regionid))
        if is_2regions_diff:
            for regionid in posinfo["region_ids"]:
                regions_diff.append((fname_prefix, regionid))
        if is_2regions_diff_sig:
            for regionid in posinfo["region_ids"]:
                regions_diff_sig.append((fname_prefix, regionid))

        if args.outdir is not None:
            wfile_all = args.outdir + "/" + fname_prefix + ".repeat_region.poses_all.txt"
            wfile_diff = args.outdir + "/" + fname_prefix + ".repeat_region.poses_diff.txt"
            if len(poses_all) > 0:
                write_poses_to_file(poses_all, wfile_all)
            if len(poses_diff) > 0:
                write_poses_to_file(poses_diff, wfile_diff)

    # write posinfo in BED format (2 kind of chrom names)
    # single posinfo / repeat info
    # (chrom, p, strand, rgid, repeat_group_id)
    fname = os.path.basename(args.region_dir.rstrip("/"))
    fname = args.outdir + "/" + fname
    # all
    allposes_bed_single_file = fname + ".posinfo.all_poses.bed"
    write_poses_to_file_bed(allrepeats_poses_all, allposes_bed_single_file)
    allposes_bed_region_file = fname + ".posinfo.all_regions.bed"
    write_posregions_to_file_bed(allrepeats_poses_all, allposes_bed_region_file)
    # all 2 group
    all2poses_bed_single_file = fname + ".posinfo.all2_poses.bed"
    write_poses_to_file_bed(allrepeats_poses_all2, all2poses_bed_single_file)
    all2poses_bed_region_file = fname + ".posinfo.all2_regions.bed"
    write_posregions_to_file_bed(allrepeats_poses_all2, all2poses_bed_region_file)
    # diff
    diffposes_bed_single_file = fname + ".posinfo.diff_poses.bed"
    write_poses_to_file_bed(allrepeats_poses_diff, diffposes_bed_single_file)

    poses_diff_in_diffregions = filter_poses_by_regionid(allrepeats_poses_diff, regions_diff)
    diffposes_indiffrgs_bed_single_file = fname + ".posinfo.diff_poses_in_diffregion.bed"
    write_poses_to_file_bed(poses_diff_in_diffregions, diffposes_indiffrgs_bed_single_file)

    poses_diff_in_sigdiffregions = filter_poses_by_regionid(allrepeats_poses_diff, regions_diff_sig)
    diffposes_insigdiffrgs_bed_single_file = fname + ".posinfo.diff_poses_in_sigdiffregion.bed"
    write_poses_to_file_bed(poses_diff_in_sigdiffregions, diffposes_insigdiffrgs_bed_single_file)

    diffposes_bed_region_file = fname + ".posinfo.diff_regions.bed"
    rg2methylevel = write_posregions_to_file_bed(allrepeats_poses_diff, diffposes_bed_region_file)

    # regions
    group_cmp_file = fname + ".group.2repeat_g1C_cov05.txt"
    write_regions_to_file(regions_cmp, group_cmp_file)

    group_diff_file = fname + ".group.2repeat_diff.txt"
    write_regions_to_file(regions_diff, group_diff_file)
    group_diff_bed = fname + ".group.2repeat_diff.bed"
    write_regions_to_bed(regions_diff, group_diff_bed, rg2methylevel)

    group_diffsig_file = fname + ".group.2repeat_diff_sig.txt"
    write_regions_to_file(regions_diff_sig, group_diffsig_file)
    group_diffsig_bed = fname + ".group.2repeat_diff_sig.bed"
    write_regions_to_bed(regions_diff_sig, group_diffsig_bed, rg2methylevel)

    # write info of repeat_groups with 2 repeats
    all2_repeats_stats_file = fname + ".all2_repeat_group_stats.txt"
    with open(all2_repeats_stats_file, "w") as wf:
        wf.write(all2_repeats_stats)

    print("repeat_group with 2 repeats, len>=100: {} ,\n"
          "    len>=1k: {},\n"
          "    len>=10k: {}".format(cnt_100, cnt_1k, cnt_10k))
    print("g\tg_2reps\tg_2reps_>=1Cs\tg_2reps_cmp\tg_2reps_diff\tg_2reps_diff_sig")
    print("{}\t{}\t{}\t{}\t{}\t{}".format(cnt_all, cnt_all2, cnt_all2_1mCs, cnt_all2_cmp,
                                          cnt_all2_diff, cnt_all2_diff_sig))
    print("=====")
    print("sites in all n_repeats_group: {} sites/{} groups;\n"
          "sites in all 2_repeats_group: {} sites/{} groups;\n"
          "diff sites in 2_repeats_group: {} sites/{} groups.".format(len(allrepeats_poses_all), cnt_all,
                                                                      len(allrepeats_poses_all2), cnt_all2,
                                                                      len(allrepeats_poses_diff), cnt_diff))


if __name__ == '__main__':
    main()
