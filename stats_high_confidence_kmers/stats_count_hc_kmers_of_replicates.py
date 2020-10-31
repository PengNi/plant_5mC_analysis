#! /usr/bin/python
import argparse

sep = "\t"


basepairs = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
             'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y',
             'Y': 'R', 'B': 'V', 'V': 'B', 'D': 'H', 'H': "D",
             'Z': 'Z'}

iupac_alphabets = {'A': ['A'], 'T': ['T'], 'C': ['C'], 'G': ['G'],
                   'R': ['A', 'G'], 'M': ['A', 'C'], 'S': ['C', 'G'],
                   'Y': ['C', 'T'], 'K': ['G', 'T'], 'W': ['A', 'T'],
                   'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
                   'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
                   'N': ['A', 'C', 'G', 'T']}


def str2bool(v):
    # susendberg's function
    return v.lower() in ("yes", "true", "t", "1")


def complement_seq(dnaseq):
    rdnaseq = dnaseq[::-1]
    comseq = ''
    try:
        comseq = ''.join([basepairs[x] for x in rdnaseq])
    except Exception:
        print('something wrong in the dna sequence.')
    return comseq


def read_position_file(positionfp, header=False):
    poses = set()
    with open(positionfp, 'r') as rf:
        if header:
            next(rf)
        for line in rf:
            words = line.strip().split("\t")
            poses.add((words[0], int(words[1]), words[2]))
    return poses


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


def _get_kmer2posinfo(contigs, poses, neighbor_len):
    kmer2posinfo = dict()
    kmerkeys = set()
    for posinfo in poses:
        chrom, pos, strand = posinfo
        pleft, pright = pos-neighbor_len, pos+neighbor_len+1
        if pleft >= 0 and pright <= len(contigs[chrom]):
            kmer = contigs[chrom][pleft:pright]
            if strand == "-":
                kmer = complement_seq(kmer)
            if kmer not in kmerkeys:
                kmer2posinfo[kmer] = set()
                kmerkeys.add(kmer)
            kmer2posinfo[kmer].add((chrom, pos, strand))
    return kmer2posinfo


def _write_positions(poses, wfp):
    poses = sorted(list(poses))
    with open(wfp, 'w') as wf:
        for pos in poses:
            wf.write("\t".join(list(map(str, pos))) + "\n")


def stats_hc_kmers(args):
    ref = DNAReference(args.ref_fa)
    contigs = ref.getcontigs()
    del ref

    hcposes1 = read_position_file(args.hcloc1, str2bool(args.header))
    hcposes2 = read_position_file(args.hcloc2, str2bool(args.header))
    print("hc poses 1: {}, hc poses 2: {}".format(len(hcposes1),
                                                  len(hcposes2)))

    klen = args.klen

    neighbor_len = klen // 2
    kmermx2poses1 = _get_kmer2posinfo(contigs, hcposes1, neighbor_len)
    kmermx2poses2 = _get_kmer2posinfo(contigs, hcposes2, neighbor_len)

    kmers1 = set(kmermx2poses1.keys())
    kmers2 = set(kmermx2poses2.keys())

    print("{}-mer, hc1: {}, hc2: {}, intersection: {}".format(klen,
                                                              len(kmers1),
                                                              len(kmers2),
                                                              len(kmers1.intersection(kmers2))))


def main():
    # python stats_count_hc_kmers_of_replicates.py --ref_fa ~/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna
    #  --hcloc1 ../stats_high_confidence_sites/ninanjie-2.bs_3replicates.CG.hc_sites.main_genome.stat.pos.inter.txt
    # --hcloc2 ../stats_high_confidence_sites/ninanjie-2.bs_3replicates.CG.hc_sites.main_genome.stat.neg.inter.txt
    # --header no --klen 13 > ninanjie-2.bs_3replicates.CG.hc_sites.main_genome.pos_inter_vs_neg_inter_kmer.log
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref_fa", type=str, required=True, help="ref path")
    parser.add_argument("--hcloc1", type=str, required=True, help="hc loc 1")
    parser.add_argument("--hcloc2", type=str, required=True, help="hc loc 2")
    parser.add_argument("--header", type=str, required=False, default="no",
                        help="pos file with header or not")
    parser.add_argument("--klen", type=int, required=False, default=13, help="kmer len")

    args = parser.parse_args()
    stats_hc_kmers(args)


if __name__ == '__main__':
    main()
