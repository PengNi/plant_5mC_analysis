#! /usr/bin/python
import argparse
import math
import os


chromname_map_arab = {"NC_003070.9": "chr1",
                      "NC_003071.7": "chr2",
                      "NC_003074.8": "chr3",
                      "NC_003075.7": "chr4",
                      "NC_003076.8": "chr5",
                      "NC_037304.1": "chrM",
                      "NC_000932.1": "chrC"}


def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")


class DNAReference:
    def __init__(self, reffile):
        self._contignames = []
        self._contigs = {}  # contigname 2 contigseq
        with open(reffile, 'r') as rf:
            contigname = ''
            contigseq = ''
            for line in rf:
                if line.startswith('>'):
                    if contigname != '' and contigseq != '':
                        self._contigs[contigname] = contigseq
                        self._contignames.append(contigname)
                    contigname = line.strip()[1:].split(' ')[0]
                    contigseq = ''
                else:
                    # turn to upper case
                    contigseq += line.strip().upper()
            self._contigs[contigname] = contigseq
            self._contignames.append(contigname)

    def getcontigs(self):
        return self._contigs

    def getcontignames(self):
        return self._contignames


def read_one_posfile(posfile):
    motif2posinfo = dict()
    with open(posfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, motif = words[0], int(words[1]), words[2], words[3]
            if motif not in motif2posinfo.keys():
                motif2posinfo[motif] = []
            motif2posinfo[motif].append((chrom, pos, strand))
    return motif2posinfo


def get_contig_range(contiglen, binlen):
    r_start = []
    for i in range(0, contiglen, binlen):
        r_start.append(i)
    r_end = r_start[1:] + [contiglen, ]
    return zip(r_start, r_end)


def get_contig_info(contigname, contigseq, binlen):
    contiglen = len(contigseq)
    contig_range = get_contig_range(contiglen, binlen)
    contigchunk2info = {}
    for cchunk in contig_range:
        r_start, r_end = cchunk
        # chuck, cg, chg, chh
        contigchunk2info[(contigname, r_start)] = [cchunk, 0, 0, 0]
    return contigchunk2info


def get_genome_contig_info(contigname2contigseq, binlen, species):
    contigchunk2info = {}
    for contigname in contigname2contigseq.keys():
        if species == "arab" and not contigname.startswith("NC_003"):
            continue
        if species == "rice" and contigname not in {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"}:
            continue
        contigchunk2info.update(get_contig_info(contigname, contigname2contigseq[contigname],
                                                binlen))
    return contigchunk2info


def fillup_contigchuninfo(contigchunk2info, motif2pos, binlen, islog):
    motif2index = {"CG": 1, "CHG": 2, "CHH": 3}
    for motif in motif2pos.keys():
        posinfo = motif2pos[motif]
        motifidx = motif2index[motif]
        for pos in posinfo:
            chrom, loc, strand = pos
            try:
                rstart_idx = loc - (loc % binlen)
                contigchunk2info[(chrom, rstart_idx)][motifidx] += 1
            except KeyError:
                print("key error")
    if islog:
        for contigchunk in contigchunk2info.keys():
            # cchunk, cg_num, chg_num, chh_num
            cchunk, cgnum, chg_num, chh_num = contigchunk2info[contigchunk]
            if cgnum != 0:
                cgnum = math.log(cgnum, 2)
            if chg_num != 0:
                chg_num = math.log(chg_num, 2)
            if chh_num != 0:
                chh_num = math.log(chh_num, 2)
            contigchunk2info[contigchunk] = cchunk, cgnum, chg_num, chh_num
    return contigchunk2info


def write_contigchunkinfo(contigchunk2info, file_prefix, species):
    methods = {1: 'cg', 2: 'chg', 3: 'chh'}
    chunkkeys = sorted(list(contigchunk2info.keys()))
    for midx in methods.keys():
        filepath = file_prefix + '.' + str(methods[midx]) + '.txt'
        nummax = 0
        with open(filepath, 'w') as wf:
            for contigchunk in chunkkeys:
                ccinfo = contigchunk2info[contigchunk]
                if nummax < ccinfo[midx]:
                    nummax = ccinfo[midx]
                chromname = chromname_map_arab[contigchunk[0]] if species == "arab" else contigchunk[0]
                wf.write(' '.join([chromname, str(ccinfo[0][0]), str(ccinfo[0][1]),
                                   str(ccinfo[midx])]) + '\n')
        print(methods[midx], nummax)


def main():
    # python circos_histogram.data.py --cytosine_info ninanjie.c_count.bs_3reps_vs_nano_guppy_part2_50x.only_nanopore_detected_Cs.pos.txt --species arab --ref_fp ~/tools/data/arab/GCF_000001735.4_TAIR10.1_genomic.fna --bin_len 100000 &
    parser = argparse.ArgumentParser("")
    parser.add_argument("--cytosine_info", type=str, required=True,
                        help="cytosine location file, from compare_detected_cytosines.py")
    parser.add_argument("--species", type=str, required=True, choices=["arab", "rice"], help="")
    parser.add_argument('-r', "--ref_fp", help="genome reference",
                        type=str, required=True)
    parser.add_argument('--bin_len', type=int, default=1000000,
                        required=False, help='bin range length')
    parser.add_argument('--islog', type=str, default='no',
                        required=False, help='is log scale')

    args = parser.parse_args()

    contigs = DNAReference(args.ref_fp).getcontigs()
    motif2posinfo = read_one_posfile(args.cytosine_info)
    chuck2info = get_genome_contig_info(contigs, args.bin_len, args.species)
    chuck2info = fillup_contigchuninfo(chuck2info, motif2posinfo, args.bin_len, str2bool(args.islog))

    fname, fext = os.path.splitext(args.cytosine_info)
    fileprefix = fname + '.circos_' + str(args.bin_len)
    if str2bool(args.islog):
        fileprefix += '.log'
    write_contigchunkinfo(chuck2info, fileprefix, args.species)


if __name__ == '__main__':
    main()
