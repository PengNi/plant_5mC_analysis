#! /usr/bin/python
import argparse

alphabeta = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def _alphabeta(letter):
    if letter in alphabeta.keys():
        return alphabeta[letter]
    return 'N'


def complement_seq(dnaseq):
    rdnaseq = dnaseq[::-1]
    comseq = ''
    try:
        comseq = ''.join([_alphabeta(x) for x in rdnaseq])
    except Exception:
        print('something wrong in the dna sequence.')
    return comseq


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


def get_refloc_of_methysite_in_motif(seqstr, motifset, methyloc_in_motif=0):
    """

    :param seqstr:
    :param motifset:
    :param methyloc_in_motif: 0-based
    :return:
    """
    motifset = set(motifset)
    strlen = len(seqstr)
    motiflen = len(list(motifset)[0])
    sites = []
    for i in range(0, strlen - motiflen + 1):
        if seqstr[i:i + motiflen] in motifset:
            sites.append(i+methyloc_in_motif)
    return sites


def convert_motif_seq(ori_seq):
    alphabets = {'A': ['A', ], 'T': ['T', ], 'C': ['C', ], 'G': ['G', ],
                 'N': ['N', ], 'H': ['A', 'C', 'T']}
    outbases = []
    for bbase in ori_seq:
        outbases.append(alphabets[bbase])

    def recursive_permute(bases_list):
        if len(bases_list) == 1:
            return bases_list[0]
        if len(bases_list) == 2:
            pseqs = []
            for fbase in bases_list[0]:
                for sbase in bases_list[1]:
                    pseqs.append(fbase + sbase)
            return pseqs
        else:
            pseqs = recursive_permute(bases_list[1:])
            pseq_list = [bases_list[0], pseqs]
            return recursive_permute(pseq_list)
    return recursive_permute(outbases)


def _count_all_motif_numbers(genomefp, motifs_str, meloc, wfile, contig_pre, contig_names):
    motifs = motifs_str.strip().split(',')
    motif_seqs = []
    for ori_motif in motifs:
        motif_seqs.append(convert_motif_seq(ori_motif))

    contigs = DNAReference(genomefp).getcontigs()
    motif_cnt_in_genome = [0 for _ in motif_seqs]
    if contig_pre is not None or contig_names is not None:
        motif_cnt_contigs = [0 for _ in motif_seqs]
    wf = open(wfile, "w")
    wf.write("contig\t" + "\t".join(motifs) + "\n")
    for contig in contigs.keys():
        ctmpseq = contigs[contig]
        ctmpseq_rc = complement_seq(ctmpseq)
        # ctmplen = len(ctmpseq)

        motifcnt = []
        # each motif
        for motif in motif_seqs:
            # forward strand
            f_poses = get_refloc_of_methysite_in_motif(ctmpseq, motif, meloc)
            # backward strand
            b_poses = get_refloc_of_methysite_in_motif(ctmpseq_rc, motif, meloc)
            print('there are {} motifs in contig {}'.format(len(f_poses) + len(b_poses), contig))
            motifcnt.append(len(f_poses) + len(b_poses))
        wf.write(contig + "\t" + "\t".join(list(map(str, motifcnt))) + "\n")
        for i in range(0, len(motif_seqs)):
            motif_cnt_in_genome[i] += motifcnt[i]
        if contig_pre is not None and str(contig).startswith(contig_pre):
            for i in range(0, len(motif_seqs)):
                motif_cnt_contigs[i] += motifcnt[i]
        if contig_names is not None:
            contigset = set(contig_names.strip().split(","))
            if contig in contigset:
                for i in range(0, len(motif_seqs)):
                    motif_cnt_contigs[i] += motifcnt[i]
    wf.write("genome" + "\t" + "\t".join(list(map(str, motif_cnt_in_genome))) + "\n")
    if contig_pre is not None or contig_names is not None:
        wf.write("main_contigs" + "\t" + "\t".join(list(map(str, motif_cnt_contigs))) + "\n")
    wf.close()


def main():
    parser = argparse.ArgumentParser("count motif num of each contig and genome")
    parser.add_argument('-r', "--ref_fp", help="the directory of genome reference",
                        type=str, required=True)
    parser.add_argument('-m', "--motifs", default='CG',
                        help="motifs, splited by comma",
                        type=str, required=False)
    parser.add_argument('--mloc_in_motif', type=int, required=False,
                        default=0,
                        help='0-based location of the methylation base in the motif, default 0')
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12
    parser.add_argument("--wfile", type=str, required=False, default=None)
    argv = parser.parse_args()
    ref_fp = argv.ref_fp
    motifs = argv.motifs
    mloc = argv.mloc_in_motif
    contig_pre = argv.contig_prefix
    contig_names = argv.contig_names
    wfile = argv.wfile
    if wfile is None:
        wfile = ref_fp + "motif_count.txt"

    _count_all_motif_numbers(ref_fp, motifs, mloc, wfile, contig_pre, contig_names)


if __name__ == '__main__':
    main()
