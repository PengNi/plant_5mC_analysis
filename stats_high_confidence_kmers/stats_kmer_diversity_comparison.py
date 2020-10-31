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


# def findsubsets(s, n):
# https://www.geeksforgeeks.org/python-program-to-get-all-subsets-of-given-size-of-a-set/
def findsubsets(s):
    import itertools
    slen = len(s)
    subsets = []
    for n in range(1, slen+1):
        subsets += [list(i) for i in itertools.combinations(s, n)]
    return subsets


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


def count_intersection(args):
    ref = DNAReference(args.ref_fa)
    contigs = ref.getcontigs()
    del ref

    klen = args.klen
    neighbor_len = klen // 2

    posfiles = args.posfile
    poses = []
    fcnt = 0
    for posfile in posfiles:
        print("{}:{}".format(fcnt, posfile))
        fcnt += 1
        poses.append(read_position_file(posfile, str2bool(args.header)))
    print("intersection stats=====")
    print("sets:pos:kmer")

    subsets = findsubsets([i for i in range(0, len(posfiles))])
    for subset in subsets:
        if len(subset) == 1:
            postmp = poses[subset[0]]
            interposes_len = len(postmp)
            kmermx2poses = _get_kmer2posinfo(contigs, postmp, neighbor_len)
            interkmers_len = len(set(kmermx2poses.keys()))
        else:
            posestmp = []
            kmerstmp = []
            for idx in subset:
                posestmp.append(poses[idx])
                kmermx2poses = _get_kmer2posinfo(contigs, poses[idx], neighbor_len)
                kmerstmp.append(set(kmermx2poses.keys()))
            posinter = set.intersection(*posestmp)
            kmerinter = set.intersection(*kmerstmp)

            interposes_len = len(posinter)
            interkmers_len = len(kmerinter)
        print("{}:{}:{}".format(subset, interposes_len, interkmers_len))


def main():
    # python stats_kmer_diversity_comparison.py
    # --ref_fa ~/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna
    # --posfile ../stats_high_confidence_sites/
    # ninanjie-2.bs_3replicates.CHG.hc_sites.main_genome.stat.pos.union.txt
    # --posfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CHG.bn11.sn128.signal_feas.bs.
    # hc_poses_all.negkmeraspos_20m.denoise6.negkmeraspos.label1.tsv.poses
    # --posfile ../stats_high_confidence_sites/
    # ninanjie-2.bs_3replicates.CHG.hc_sites.main_genome.stat.pos.inter.txt
    # > ninanjie-2.bs_3replicates.CHG.hc_sites.main_genome.pos_kmer_diversity.log
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref_fa", type=str, required=True, help="ref path")
    parser.add_argument("--posfile", type=str, action="append", required=True, help="")
    parser.add_argument("--header", type=str, required=False, default="no",
                        help="pos file with header or not")
    parser.add_argument("--klen", type=int, required=False, default=11, help="kmer len")

    args = parser.parse_args()
    count_intersection(args)


if __name__ == '__main__':
    main()
