#! /usr/bin/python
import argparse
import os

sep = "||"
chromname_map_arab = {"NC_003070.9": "Chr1",
                      "NC_003071.7": "Chr2",
                      "NC_003074.8": "Chr3",
                      "NC_003075.7": "Chr4",
                      "NC_003076.8": "Chr5",
                      "NC_037304.1": "ChrM",
                      "NC_000932.1": "ChrC"}
chromname_map_rice = {"NC_029256.1": "1",
                      "NC_029257.1": "2",
                      "NC_029258.1": "3",
                      "NC_029259.1": "4",
                      "NC_029260.1": "5",
                      "NC_029261.1": "6",
                      "NC_029262.1": "7",
                      "NC_029263.1": "8",
                      "NC_029264.1": "9",
                      "NC_029265.1": "10",
                      "NC_029266.1": "11",
                      "NC_029267.1": "12"}


# ======================
def read_cytosine_info(cytosine_info_file):
    posinfo = []
    with open(cytosine_info_file, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            chrom, loc, strand, motif = words[0], int(words[1]), words[2], words[-1]
            posinfo.append((chrom, loc, strand, motif))
    return posinfo


def stat_posinfo(posinfo):
    motif2num = {}
    for tpos in posinfo:
        chrom, loc, strand, motif = tpos
        if motif not in motif2num.keys():
            motif2num[motif] = 0
        motif2num[motif] += 1
    return motif2num


def reformat_posinfo_for_cmp(posinfo):
    chrom2motifsites = {}
    for tposinfo in posinfo:
        chrom, loc, strand, motif = tposinfo
        # keytmp = sep.join([chrom, strand])
        keytmp = chrom
        if keytmp not in chrom2motifsites.keys():
            chrom2motifsites[keytmp] = {}
        if motif not in chrom2motifsites[keytmp].keys():
            chrom2motifsites[keytmp][motif] = []
        chrom2motifsites[keytmp][motif].append(loc)
    return chrom2motifsites


# read trf/irf .dat file
def read_regioninfo_from_datfile(datfile):
    regions = []
    with open(datfile, 'r') as rf:
        chromname = ""
        for line in rf:
            if line.startswith("Sequence:"):
                chromname = line.strip().split(" ")[1]
            elif line[0].isdigit():
                words = line.strip().split(" ")
                start, end = int(words[0]) - 1, int(words[1])
                regions.append((chromname, start, end))
    return regions


def read_regioninfo_from_bedfile(bedfile):
    regions = []
    with open(bedfile, 'r') as rf:
        for line in rf:
            if not line.startswith("#"):
                words = line.strip().split("\t")
                chrom, start, end = words[0], int(words[1]), int(words[2])
                regions.append((chrom, start, end))
    return regions


def reformat_regioninfo_to_chrom2regions(regioninfo, species="arab", is_from_ncbi=False):
    chrom2regions = dict()
    for region in regioninfo:
        chrom, start, end = region
        if species == "rice" and is_from_ncbi:
            chrom = chromname_map_rice[chrom]
        if chrom not in chrom2regions.keys():
            chrom2regions[chrom] = []
        chrom2regions[chrom].append((start, end))
    return chrom2regions


def flat_region_sites(regions):
    regionsites = set()
    for region in regions:
        start, end = region
        for locidx in range(start, end):
            regionsites.add(locidx)
    return regionsites


def get_pos_distribution_in_generegion(rgtype2regions, chrom2motifsites):
    rgtype2stats = {}
    for rgtype in rgtype2regions.keys():
        motif2count = {}
        rgsitesnum = 0
        chrom2regions = rgtype2regions[rgtype]

        for chromkey in chrom2regions.keys():
            if chromkey in chrom2motifsites.keys():
                # print(rgtype, chromkey)
                regionsites = set()
                for region in chrom2regions[chromkey]:
                    start, end = region
                    for locidx in range(start, end):
                        regionsites.add(locidx)
                rgsitesnum += len(regionsites)

                motif2sites = chrom2motifsites[chromkey]
                for motif in motif2sites.keys():
                    sitestmp = set(motif2sites[motif])
                    num_sites_in_region = len(sitestmp.intersection(regionsites))
                    if motif not in motif2count.keys():
                        motif2count[motif] = 0
                    motif2count[motif] += num_sites_in_region
        print("{} sitesnum: {}, motif-intersect-{}".format(rgtype, rgsitesnum, sum(list(motif2count.values()))))
        rgtype2stats[rgtype] = motif2count
    return rgtype2stats


def stat_cytosine_distribution_in_repeats(args):
    fname, fext = os.path.splitext(args.cytosine_info)
    posinfo = read_cytosine_info(args.cytosine_info)
    motif2counts = stat_posinfo(posinfo)
    print(motif2counts)
    cytosines_nano = reformat_posinfo_for_cmp(posinfo)
    totalnum = len(posinfo)

    repeat_window = read_regioninfo_from_bedfile(args.windowmasker_bed)
    print("windowmasker region num: {}".format(len(repeat_window)))
    regions_win = reformat_regioninfo_to_chrom2regions(repeat_window, args.species, True)
    repeat_repeat = read_regioninfo_from_bedfile(args.repeatmasker_bed)
    print("repeatmasker region num: {}".format(len(repeat_repeat)))
    regions_rep = reformat_regioninfo_to_chrom2regions(repeat_repeat, args.species, True)
    repeat_trf = read_regioninfo_from_datfile(args.trf_dat)
    regions_trf = reformat_regioninfo_to_chrom2regions(repeat_trf, args.species)
    repeat_irf = read_regioninfo_from_datfile(args.irf_dat)
    regions_irf = reformat_regioninfo_to_chrom2regions(repeat_irf, args.species)

    rg2regions = dict()
    rg2regions["Repeats(WindowMasker)"] = regions_win
    rg2regions["Repeats(RepeatMasker)"] = regions_rep
    rg2regions["Tandem repeats"] = regions_trf
    rg2regions["Inverted repeats"] = regions_irf
    # stat
    rg2stats = get_pos_distribution_in_generegion(rg2regions, cytosines_nano)
    wfile = fname + ".cnum_in_repeats.txt"
    with open(wfile, "w") as wf:
        wf.write("\t".join(["region", "motif", "moitf_num", "motif_total", "c_total"]) + "\n")
        for rgtype in rg2stats.keys():
            for motif in rg2stats[rgtype].keys():
                wf.write("\t".join([rgtype, motif, str(rg2stats[rgtype][motif]),
                                    str(motif2counts[motif]), str(totalnum)]) + "\n")


def main():
    # python generate_Cdistri_in_repeats.py --cytosine_info ninanjie.c_count.bs_3reps_vs_nano50x.only_nanopore_detected_Cs.pos.txt
    # --windowmasker_bed data.arab/Repeats_identified_by_WindowMasker.chr1-5.BED
    # --trf_dat data.arab/trf.GCF_000001735.4_TAIR10.1_genomic.fna.2.7.7.80.10.50.500.dat
    # --irf_dat data.arab/irf.GCF_000001735.4_TAIR10.1_genomic.fna.2.3.5.80.10.40.500000.10000.dat
    # --repeatmasker_bed data.arab/Repeats_identified_by_RepeatMasker.chr1-5.BED &
    parser = argparse.ArgumentParser("repeat regions")
    parser.add_argument("--cytosine_info", type=str, required=True,
                        help="cytosine location file, from generate_nanopore_only_cytosines_info.py")
    parser.add_argument("--species", type=str, required=False, default="arab", choices=["arab", "rice"], help="")
    parser.add_argument("--windowmasker_bed", type=str, required=True,
                        help="")
    parser.add_argument("--repeatmasker_bed", type=str, required=True,
                        help="")
    parser.add_argument("--trf_dat", type=str, required=True,
                        help="")
    parser.add_argument("--irf_dat", type=str, required=True,
                        help="")

    args = parser.parse_args()
    stat_cytosine_distribution_in_repeats(args)


def test():
    repeat_window = read_regioninfo_from_bedfile("data.arab/Repeats_identified_by_WindowMasker.chr1-5.BED")
    repeat_repeat = read_regioninfo_from_bedfile("data.arab/Repeats_identified_by_RepeatMasker.chr1-5.BED")
    regions_win = reformat_regioninfo_to_chrom2regions(repeat_window)
    regions_rep = reformat_regioninfo_to_chrom2regions(repeat_repeat)

    for chrom in regions_win.keys():
        if chrom in regions_rep.keys():
            win_sites = flat_region_sites(regions_win[chrom])
            rep_sites = flat_region_sites(regions_rep[chrom])
            print("chrom-{}: windowmasker-{}, repeatmasker-{}, inersect-{}".format(chrom,
                                                                                   len(win_sites),
                                                                                   len(rep_sites),
                                                                                   len(win_sites.intersection(rep_sites))))

    print("======")
    repeat_trf = read_regioninfo_from_datfile("data.arab/trf.GCF_000001735.4_TAIR10.1_genomic.fna.2.7.7.80.10.50.500.dat")
    regions_trf = reformat_regioninfo_to_chrom2regions(repeat_trf)
    print("trf len: {}".format(len(repeat_trf)))
    for chrom in regions_trf.keys():
        if chrom in regions_rep.keys():
            trf_sites = flat_region_sites(regions_trf[chrom])
            win_sites = flat_region_sites(regions_win[chrom])
            rep_sites = flat_region_sites(regions_rep[chrom])
            print("chrom-{}: trf-{}, intersect_repeatmasker-{}, inersect_windowmasker-{}".format(chrom,
                                                                                                 len(trf_sites),
                                                                                                 len(trf_sites.intersection(rep_sites)),
                                                                                                 len(trf_sites.intersection(win_sites))))

    print("======")
    repeat_irf = read_regioninfo_from_datfile(
        "data.arab/irf.GCF_000001735.4_TAIR10.1_genomic.fna.2.3.5.80.10.40.500000.10000.dat")
    regions_irf = reformat_regioninfo_to_chrom2regions(repeat_irf)
    print("irf len: {}".format(len(repeat_irf)))
    for chrom in regions_trf.keys():
        if chrom in regions_rep.keys():
            irf_sites = flat_region_sites(regions_irf[chrom])
            win_sites = flat_region_sites(regions_win[chrom])
            rep_sites = flat_region_sites(regions_rep[chrom])
            print("chrom-{}: irf-{}, intersect_repeatmasker-{}, inersect_windowmasker-{}".format(chrom,
                                                                                                 len(irf_sites),
                                                                                                 len(
                                                                                                     irf_sites.intersection(
                                                                                                         rep_sites)),
                                                                                                 len(
                                                                                                     irf_sites.intersection(
                                                                                                         win_sites))))


if __name__ == '__main__':
    main()
    # test()
