#! /usr/bin/python
import argparse

key_sep = "||"


def _read_one_bs_rmet_file(bs_file, args):
    hrmetsets = dict()
    lrmetsets = dict()
    irmetsets = dict()
    contigset = set(args.contig_names.strip().split(",")) if args.contig_names is not None else None
    with open(bs_file, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, coverage, rmet = words[0], words[1], words[2], int(words[6]), float(words[7])
            mkey = key_sep.join([chrom, str(pos)])
            cnt_flag = 0
            if args.contig_prefix is not None:
                if str(chrom).startswith(args.contig_prefix):
                    cnt_flag = 1
            elif args.contig_names is not None:
                if chrom in contigset:
                    cnt_flag = 1
            else:
                cnt_flag = 1
            if cnt_flag == 1:
                if coverage >= args.covcf:
                    if rmet >= args.hrmet:
                        hrmetsets[mkey] = rmet
                    elif rmet <= args.lrmet:
                        lrmetsets[mkey] = rmet
                    else:
                        irmetsets[mkey] = rmet
    print("{}: hrmet,{}; lrmet,{}".format(bs_file, len(hrmetsets), len(lrmetsets)))
    return hrmetsets, lrmetsets, irmetsets


def _read_bs_rmet_files(args):
    print("==bs rmet files==:")
    bsfiles = args.bsfile
    bsfiles_len = len(bsfiles)
    for idx in range(0, bsfiles_len):
        print("=={}, {}".format(idx, bsfiles[idx]))
    print("covcf: {}, hrmet: {}, lrmet: {}, main_genome: {}/{}".format(args.covcf,
                                                                       args.hrmet,
                                                                       args.lrmet,
                                                                       args.contig_names,
                                                                       args.contig_prefix))

    hrmetkeys, lrmetkeys, irmetkeys = [], [], []
    hrmetsets, lrmetsets, irmetsets = [], [], []
    for bsfile in bsfiles:
        hrmetset, lrmetset, irmetset = _read_one_bs_rmet_file(bsfile, args)
        hrmetsets.append(hrmetset)
        lrmetsets.append(lrmetset)
        irmetsets.append(irmetset)
        hrmetkeys.append(set(hrmetset.keys()))
        lrmetkeys.append(set(lrmetset.keys()))
        irmetkeys.append(set(irmetset.keys()))
    hrmetkeys = set.intersection(*hrmetkeys)
    lrmetkeys = set.intersection(*lrmetkeys)
    irmetkeys = set.intersection(*irmetkeys)
    print("==bs replicates intersection: hrmet,{}; lrmet,{}; irmet,{}".format(len(hrmetkeys), len(lrmetkeys),
                                                                              len(irmetkeys)))
    avg_hrmet, avg_lrmet, avg_irmet = dict(), dict(), dict()
    for hrmetkey in hrmetkeys:
        rmets = []
        for hrmetset in hrmetsets:
            rmets.append(hrmetset[hrmetkey])
        rmet_avg = sum(rmets) / float(len(rmets))
        avg_hrmet[hrmetkey] = rmet_avg
    for lrmetkey in lrmetkeys:
        rmets = []
        for lrmetset in lrmetsets:
            rmets.append(lrmetset[lrmetkey])
        rmet_avg = sum(rmets) / float(len(rmets))
        avg_lrmet[lrmetkey] = rmet_avg
    for irmetkey in irmetkeys:
        rmets = []
        for irmetset in irmetsets:
            rmets.append(irmetset[irmetkey])
        rmet_avg = sum(rmets) / float(len(rmets))
        avg_irmet[irmetkey] = rmet_avg
    return avg_hrmet, avg_lrmet, avg_irmet


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bsfile", action="append", type=str, required=True,
                        help="bs rmet file")
    parser.add_argument("--wfile", type=str, required=True,
                        help="")
    parser.add_argument("--covcf", type=int, required=False, default=5,
                        help="")
    parser.add_argument("--hrmet", type=float, required=False, default=0.7,
                        help="")
    parser.add_argument("--lrmet", type=float, required=False, default=0.3,
                        help="")
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12

    args = parser.parse_args()

    avg_hrmet, avg_lrmet, avg_irmet = _read_bs_rmet_files(args)
    wf = open(args.wfile, 'w')
    for rmetkey in avg_hrmet.keys():
        wf.write("\t".join([rmetkey, str(avg_hrmet[rmetkey]), "high", "bisulfite"]) + "\n")
    for rmetkey in avg_lrmet.keys():
        wf.write("\t".join([rmetkey, str(avg_lrmet[rmetkey]), "low", "bisulfite"]) + "\n")
    for rmetkey in avg_irmet.keys():
        wf.write("\t".join([rmetkey, str(avg_irmet[rmetkey]), "mediate", "bisulfite"]) + "\n")
    wf.close()


if __name__ == '__main__':
    main()
