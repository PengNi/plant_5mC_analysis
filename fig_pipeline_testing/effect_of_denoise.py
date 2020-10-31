#! /usr/bin/python
import argparse
from copy import deepcopy


def read_kmer2signalinfo_from_feature_file(feafile, is_negative):
    kmer2signalinfo = dict()
    with open(feafile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            chrom = words[0]
            pos = words[1]
            rid = words[4]
            sampleid = "||".join([chrom, pos, rid])
            kmer = words[6]
            smeans = [float(x) for x in words[7].split(",")]
            label = int(words[11])
            if label == 1 or (label == 0 and is_negative):
                if kmer not in kmer2signalinfo.keys():
                    kmer2signalinfo[kmer] = dict()
                if label not in kmer2signalinfo[kmer].keys():
                    kmer2signalinfo[kmer][label] = []
                kmer2signalinfo[kmer][label].append((sampleid, smeans))
    return kmer2signalinfo


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--balance_file", type=str, required=True, help="")
    parser.add_argument("--denoise_file", type=str, required=True, help="")
    parser.add_argument("--w_file", type=str, required=True, help="")

    args = parser.parse_args()

    bkmerinfo = read_kmer2signalinfo_from_feature_file(args.balance_file, True)
    dkmerinfo = read_kmer2signalinfo_from_feature_file(args.denoise_file, False)
    for bkmer in bkmerinfo.keys():
        if 1 in bkmerinfo[bkmer].keys():
            if bkmer not in dkmerinfo.keys():
                bkmerinfo[bkmer][-1] = bkmerinfo[bkmer][1]
                del bkmerinfo[bkmer][1]
            else:
                ssinfos = dkmerinfo[bkmer][1]
                sids = set()
                for ssinfo in ssinfos:
                    sids.add(ssinfo[0])
                bkmer1ssinfo = deepcopy(bkmerinfo[bkmer][1])
                for ssinfo in bkmer1ssinfo:
                    sid, sinfo = ssinfo
                    if sid not in sids:
                        if -1 not in bkmerinfo[bkmer].keys():
                            bkmerinfo[bkmer][-1] = []
                        bkmerinfo[bkmer][-1].append((sid, sinfo))
                        bkmerinfo[bkmer][1].remove(ssinfo)
    del dkmerinfo
    wf = open(args.w_file, "w")
    wf.write("\t".join(["kmer", "label", "sid", ] + ["smean" + str(x) for x in range(0, 13)]) + "\n")
    for bkmer in bkmerinfo.keys():
        if len(bkmerinfo[bkmer]) >= 3:
            print(bkmer, len(bkmerinfo[bkmer][1]), len(bkmerinfo[bkmer][-1]), len(bkmerinfo[bkmer][0]))
            if len(bkmerinfo[bkmer][0]) >= 10 and len(bkmerinfo[bkmer][1]) >= 10 and len(bkmerinfo[bkmer][-1]) >= 10:
                for idx in [1, -1, 0]:
                    for ssinfo in bkmerinfo[bkmer][idx]:
                        sid, sinfo = ssinfo
                        wf.write("\t".join([bkmer, str(idx), sid, ] + list(map(str, sinfo))) + "\n")
    wf.flush()
    wf.close()


if __name__ == '__main__':
    main()
