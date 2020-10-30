#! /usr/bin/python
import argparse


def _read_one_bs_rmet_file(bs_file, args):
    hrmetsets = dict()
    lrmetsets = dict()
    contigset = set(args.contig_names.strip().split(",")) if args.contig_names is not None else None
    with open(bs_file, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, coverage, rmet = words[0], words[1], words[2], int(words[6]), float(words[7])
            mkey = "\t".join([chrom, str(pos), strand])
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
    print("{}: hrmet,{}; lrmet,{}".format(bs_file, len(hrmetsets), len(lrmetsets)))
    return hrmetsets, lrmetsets


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

    hrmetkeys, lrmetkeys = [], []
    hrmetsets, lrmetsets = [], []
    for bsfile in bsfiles:
        hrmetset, lrmetset = _read_one_bs_rmet_file(bsfile, args)
        hrmetsets.append(hrmetset)
        lrmetsets.append(lrmetset)
        hrmetkeys.append(set(hrmetset.keys()))
        lrmetkeys.append(set(lrmetset.keys()))
    hrmetkeys = set.intersection(*hrmetkeys)
    lrmetkeys = set.intersection(*lrmetkeys)
    print("==bs replicates intersection: hrmet,{}; lrmet,{}".format(len(hrmetkeys), len(lrmetkeys)))
    avg_hrmet, avg_lrmet = dict(), dict()
    for hrmetkey in hrmetkeys:
        rmets = []
        for hrmetset in hrmetsets:
            rmets.append(hrmetset[hrmetkey])
        rmet_avg = sum(rmets) / len(rmets)
        avg_hrmet[hrmetkey] = rmet_avg
    for lrmetkey in lrmetkeys:
        rmets = []
        for lrmetset in lrmetsets:
            rmets.append(lrmetset[lrmetkey])
        rmet_avg = sum(rmets) / len(rmets)
        avg_lrmet[lrmetkey] = rmet_avg
    return avg_hrmet, avg_lrmet


def _read_one_mod_freq_file(freqfile, hrmetsets, lrmetsets, args):
    w_strs = []
    cnt_h, cnt_l = 0, 0
    cnt_nh, cnt_nl = 0, 0
    nano_hrmet, nano_lrmet = set(), set()

    contigset = set(args.contig_names.strip().split(",")) if args.contig_names is not None else None
    with open(freqfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            chrom = words[0]
            m_key = "\t".join([words[0], words[1], words[2]])
            cov = int(words[8])
            rmet = float(words[9])

            cnt_flag = 0
            if args.contig_prefix is not None:
                if str(chrom).startswith(args.contig_prefix):
                    cnt_flag = 1
            elif args.contig_names is not None:
                if chrom in contigset:
                    cnt_flag = 1
            else:
                cnt_flag = 1

            if cnt_flag == 1 and cov >= args.covcf:
                if rmet >= args.hrmet:
                    nano_hrmet.add(m_key)
                elif rmet <= args.lrmet:
                    nano_lrmet.add(m_key)
                if m_key in hrmetsets:
                    w_strs.append("\t".join([str(rmet), "highly_nano"]))
                    cnt_h += 1
                    if rmet >= args.hrmet:
                        cnt_nh += 1
                elif m_key in lrmetsets:
                    w_strs.append("\t".join([str(rmet), "lowly_nano"]))
                    cnt_l += 1
                    if rmet <= args.lrmet:
                        cnt_nl += 1
    print("=={}:\nnanopore highly,{}; nanopore lowly,{}".format(freqfile,
                                                                len(nano_hrmet),
                                                                len(nano_lrmet)))
    print("intersection with BS: highly, {}; lowly,{}".format(len(nano_hrmet.intersection(set(hrmetsets.keys()))),
                                                              len(nano_lrmet.intersection(set(lrmetsets.keys())))))
    print("difference with BS: highly, {}; lowly,{}".format(len(nano_hrmet.difference(set(hrmetsets.keys()))),
                                                            len(nano_lrmet.difference(set(lrmetsets.keys())))))
    print("BS difference with nano: highly, {}; lowly,{}".format(len(set(hrmetsets.keys()).difference(nano_hrmet)),
                                                                 len(set(lrmetsets.keys()).difference(nano_lrmet))))
    print("\n")
    print("=={}: \nbs highly sites detected by nano,{}; bs lowly detected sites by nano, {}".format(freqfile,
                                                                                                    cnt_h,
                                                                                                    cnt_l))
    print("highly in nano&bs/highly in bs: {}, {}; "
          "lowly in nano&bs/lowly in bs: {}, {}".format(cnt_nh,
                                                        float(cnt_nh)/len(hrmetsets.keys()),
                                                        cnt_nl,
                                                        float(cnt_nl)/len(lrmetsets.keys())))
    for hrmetbskey in hrmetsets.keys():
        w_strs.append("\t".join([str(hrmetsets[hrmetbskey]), "highly_bs"]))
    for lrmetbskey in lrmetsets.keys():
        w_strs.append("\t".join([str(lrmetsets[lrmetbskey]), "lowly_bs"]))
    return w_strs


def write_rmet_strs(w_strs, wfile):
    wf = open(wfile, 'w')
    for wstr in w_strs:
        wf.write(wstr + "\n")
    wf.close()


def main():
    # python comparison_highly_lowly_methylated_sites.py --contig_prefix NC_003 --bsfile /homeb/nipeng/data/nanopore/arabidopsis/20190421-NPL0867-P1-E11-H11-barcode/bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --bsfile /homeb/nipeng/data/nanopore/arabidopsis/20190421-NPL0867-P1-E11-H11-barcode/bs.poses/D1902827A-ZJ_364364.L1L2_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --bsfile /homeb/nipeng/data/nanopore/arabidopsis/20190421-NPL0867-P1-E11-H11-barcode/bs.poses/D1902828A-ZJ_363363.L2L4_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --freqfile /homeb/nipeng/data/nanopore/arabidopsis/20190421-NPL0867-P1-E11-H11-barcode/athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.50x.12345_call_mods.freq.tsv --wfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.50x.12345_call_mods.freq.high_low_in_bs.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.50x.12345_call_mods.freq.high_low_in_bs.log
    # python comparison_highly_lowly_methylated_sites.py --contig_prefix NC_003 --bsfile ~/tools/data/arab/bs.result/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --bsfile ~/tools/data/arab/bs.result/D1902827A-ZJ_364364.L1L2_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --bsfile ~/tools/data/arab/bs.result/D1902828A-ZJ_363363.L2L4_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --freqfile ~/tools/data/arab/athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.tsv --wfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.high_low_in_bs.tsv > athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.high_low_in_bs.log &
    parser = argparse.ArgumentParser()
    parser.add_argument("--bsfile", action="append", type=str, required=True,
                        help="bs rmet file")
    parser.add_argument("--freqfile", type=str, required=True,
                        help="nanopore freq file")
    parser.add_argument("--wfile", type=str, required=True,
                        help="")
    parser.add_argument("--hrmet", type=float, required=False, default=0.7,
                        help="")
    parser.add_argument("--lrmet", type=float, required=False, default=0.3,
                        help="")
    parser.add_argument("--covcf", type=int, required=False, default=5,
                        help="")
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12

    args = parser.parse_args()

    avg_hrmet, avg_lrmet = _read_bs_rmet_files(args)
    wstrs = _read_one_mod_freq_file(args.freqfile, avg_hrmet, avg_lrmet, args)
    write_rmet_strs(wstrs, args.wfile)


if __name__ == '__main__':
    main()
