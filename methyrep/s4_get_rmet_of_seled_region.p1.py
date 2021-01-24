#!/usr/bin/python
import argparse
import os
import re
import uuid


def str2bool(v):
    # susendberg's function
    return v.lower() in ("yes", "true", "t", "1")


# DNA genome ===============================================================
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


def complement_seq(dnaseq):
    rdnaseq = dnaseq[::-1]
    comseq = ''
    try:
        comseq = ''.join([basepairs[x] for x in rdnaseq])
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


def _convert_motif_seq(ori_seq, is_dna=True):
    outbases = []
    for bbase in ori_seq:
        outbases.append(iupac_alphabets[bbase])

    def recursive_permute(bases_list):
        if len(bases_list) == 1:
            return bases_list[0]
        elif len(bases_list) == 2:
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


def get_motif_seqs(motifs, is_dna=True):
    ori_motif_seqs = motifs.strip().split(',')

    motif_seqs = []
    for ori_motif in ori_motif_seqs:
        motif_seqs += _convert_motif_seq(ori_motif.strip().upper(), is_dna)
    return motif_seqs


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


# PAF alignment ===============================================
def _get_contig_len(ref_fa):
    with open(ref_fa, "r") as rf:
        refname = next(rf).strip()[1:]
        refseq = ""
        for line in rf:
            refseq += line.strip()
        return len(refseq)


def _align_minimap2_asm(ref_fa, query_fa):
    fname, fext = os.path.splitext(ref_fa)
    paf_fp = fname + ".paf"

    if _get_contig_len(ref_fa) > 500:
        cmd = "minimap2 -x asm5 -c -t 20 {ref_fa} {query_fa} > {paf}".format(ref_fa=ref_fa,
                                                                             query_fa=query_fa,
                                                                             paf=paf_fp)
    else:
        cmd = "minimap2 -x sr -c -t 20 {ref_fa} {query_fa} > {paf}".format(ref_fa=ref_fa,
                                                                           query_fa=query_fa,
                                                                           paf=paf_fp)
    os.system(cmd)
    return paf_fp


class PafRecord2:
    def __init__(self, fields):
        self._queryname = fields[0]
        self._querylen = int(fields[1])
        self._querystart = int(fields[2])
        self._queryend = int(fields[3])
        self._strand = fields[4]
        self._targetname = fields[5]
        self._targetlen = int(fields[6])
        self._targetstart = int(fields[7])
        self._targetend = int(fields[8])
        self._match = int(fields[9])
        self._alignlen = int(fields[10])
        self._mapq = int(fields[11])
        self._cigar = fields[-1].split(":")[-1]

    def identity(self):
        return float(self._match)/self._alignlen

    def alignlen(self):
        return self._alignlen

    def is_query_equal_target(self):
        if self._queryname == self._targetname and self._strand == "+" \
                and (abs(self._querystart - self._targetstart) < 1000 or
                     abs(self._queryend - self._targetend) < 1000):
            return True
        return False


def _get_pairwise_alignment_from_cigar(queryseq, targetseq, cigarseq):
    qextseq, textseq = "", ""
    pattern = re.compile(r'((\d)+(M|S|H|X|=|I|D|N))')
    it = pattern.findall(cigarseq)

    cidx_q, cidx_t = 0, 0
    for match in it:
        num = int(match[0][:-1])
        if match[0].endswith('S') or match[0].endswith('H'):
            cidx_q += num
        elif match[0].endswith('X') or match[0].endswith('=') or match[0].endswith('M'):
            qextseq += queryseq[cidx_q:(cidx_q + num)]
            textseq += targetseq[cidx_t:(cidx_t + num)]
            cidx_q += num
            cidx_t += num
        elif match[0].endswith('I'):
            qextseq += queryseq[cidx_q:(cidx_q + num)]
            textseq += "-" * num
            cidx_q += num
        elif match[0].endswith('D'):
            qextseq += "-" * num
            textseq += targetseq[cidx_t:(cidx_t + num)]
            cidx_t += num
        elif match[0].endswith('N'):
            pass
    return qextseq, textseq


def _get_match_pair_loc(query_extseq, target_extseq):
    assert(len(query_extseq) == len(target_extseq))

    qt_match_pairs = []
    cnt_m, cnt_unm, cnt_indel = 0, 0, 0
    qidx, tidx = 0, 0
    for i in range(0, len(query_extseq)):
        if query_extseq[i] == "-":
            tidx += 1
            cnt_indel += 1
        elif target_extseq[i] == "-":
            qidx += 1
            cnt_indel += 1
        else:
            if query_extseq[i] != target_extseq[i]:
                print(i, query_extseq[i], target_extseq[i])
                cnt_unm += 1
            else:
                cnt_m += 1
            qt_match_pairs.append((qidx, tidx))
            qidx += 1
            tidx += 1
    print("total len: {}, match pair len: {}(m:{}/unm:{}), indel len: {}".format(len(query_extseq),
                                                                                 len(qt_match_pairs),
                                                                                 cnt_m, cnt_unm,
                                                                                 cnt_indel))
    return qt_match_pairs


def _get_refseq(chrom, strand, start, end, contigname2seq, is_target):
    chrom_eles = chrom.split("_")
    chrom_c, chrom_s, chrom_e = "_".join(chrom_eles[:-2]), int(chrom_eles[-2]), int(chrom_eles[-1])
    subseq = contigname2seq[chrom_c][(chrom_s-1):chrom_e]
    # TODO: make sure fanzhuan and qie de shunxu, is this right
    if is_target:
        if strand == "-":
            subseq = complement_seq(subseq)
        subseq = subseq[start:end]
    else:
        # for query seq, first qie then fanzhuan is right
        subseq = subseq[start:end]
        if strand == "-":
            subseq = complement_seq(subseq)
    return subseq


def _get_absolute_loc(relative_locs, abs_start, abs_end, strand, is_target):
    if is_target:
        return [x + abs_start for x in relative_locs]
    else:
        if strand == "+":
            return [x + abs_start for x in relative_locs]
        else:
            return [abs_end - 1 - x for x in relative_locs]


def _get_absolute_loc2(relative_locs, chrom, is_target=False, targetstrand=None):
    chrom_eles = chrom.split("_")
    chrom_c, chrom_s, chrom_e = "_".join(chrom_eles[:-2]), int(chrom_eles[-2]), int(chrom_eles[-1])
    if is_target:
        if targetstrand == "+":
            return [(chrom_s - 1) + x for x in relative_locs]
        else:
            return [(chrom_e - 1) - x for x in relative_locs]
    else:
        return [(chrom_s-1) + x for x in relative_locs]


def _get_5mer(chromloc_key, siteloc, contigs):
    chrom_eles = chromloc_key.split("_")
    chrom_c, chrom_s, chrom_e, strand = "_".join(chrom_eles[:-3]), int(chrom_eles[-3]), \
                                        int(chrom_eles[-2]), chrom_eles[-1]
    motif_s = siteloc - 2
    motif_e = siteloc + 2 + 1
    if motif_s < 0 or motif_e > len(contigs[chrom_c]):
        return None
    motif_seq = contigs[chrom_c][motif_s:motif_e]
    if strand == "-":
        motif_seq = complement_seq(motif_seq)
    return motif_seq


def parse_paf(paffile, contigname2seq, ref_strand, motifs="CG", mod_loc=0):  # parse_paf and keep only motif locs
    # WARN: only for alignment between (sub)seqs from one same genome/cotig file
    target_strand = ref_strand
    target_name = None
    motiflocs = None
    querys_locs = dict()
    with open(paffile, 'r') as rf:
        for line in rf:
            words = line.strip().split("\t")
            pafrecord2 = PafRecord2(words)
            target_name = pafrecord2._targetname
            if motiflocs is None:
                target_eles = target_name.split("_")
                target_c, target_s, target_e = "_".join(target_eles[:-2]), int(target_eles[-2]), int(target_eles[-1])
                targetseq = contigname2seq[target_c][(target_s-1):target_e]
                if target_strand == "-":
                    targetseq = complement_seq(targetseq)
                motiflocs = get_refloc_of_methysite_in_motif(targetseq, set(get_motif_seqs(motifs)), mod_loc)
                if target_strand == "+":
                    motiflocs = [x + target_s - 1 for x in motiflocs]
                else:
                    motiflocs = [target_e - 1 - x for x in motiflocs]
                motiflocs = set(motiflocs)
            queryseq = _get_refseq(pafrecord2._queryname, pafrecord2._strand,
                                   pafrecord2._querystart, pafrecord2._queryend,
                                   contigname2seq, False)
            targetseq = _get_refseq(pafrecord2._targetname, target_strand,
                                    pafrecord2._targetstart, pafrecord2._targetend,
                                    contigname2seq, True)
            # print("queryseq:\t{}".format(queryseq))
            # print("targetseq:\t{}".format(targetseq))
            query_extseq, target_extseq = _get_pairwise_alignment_from_cigar(queryseq, targetseq,
                                                                             pafrecord2._cigar)
            print("match_query:\t{}".format(query_extseq))
            print("match_target:\t{}".format(target_extseq))
            qt_match_pairs = list(zip(*_get_match_pair_loc(query_extseq, target_extseq)))
            qmatchlocs, tmatchlocs = qt_match_pairs[0], qt_match_pairs[1]
            # get absolute loc in the seq for alignment
            q_abslocs, tabslocs = _get_absolute_loc(qmatchlocs, pafrecord2._querystart,
                                                    pafrecord2._queryend, pafrecord2._strand, False), \
                                  _get_absolute_loc(tmatchlocs, pafrecord2._targetstart,
                                                    pafrecord2._targetend, target_strand, True)
            # get absolute loc in the genome reference
            q_abslocs, tabslocs = _get_absolute_loc2(q_abslocs, pafrecord2._queryname, False), \
                                  _get_absolute_loc2(tabslocs, pafrecord2._targetname, True, target_strand)

            # keep only motif locs
            # print(motiflocs)
            motiflocs_tmp = motiflocs.intersection(set(tabslocs))
            tabsloc2idx = dict()
            for i in range(0, len(tabslocs)):
                tabsloc2idx[tabslocs[i]] = i
            t2q_locs = dict()
            for motifloc in motiflocs_tmp:
                t2q_locs[motifloc] = q_abslocs[tabsloc2idx[motifloc]]
            querys_locs["_".join([pafrecord2._queryname, pafrecord2._strand])] = t2q_locs

    if len(querys_locs) == 0:
        return ""
    tlocs = []
    for query_name in querys_locs.keys():
        tlocs.append(set(querys_locs[query_name].keys()))
    tlocs_inter = sorted(list(set.intersection(*tlocs)))
    # filter tlocs by seq
    tlocs_inter_tmp = []
    for i in range(0, len(tlocs_inter)):
        tsite_motif = _get_5mer("_".join([target_name, target_strand]), tlocs_inter[i], contigname2seq)
        if tsite_motif is None:
            continue
        is_keep = True
        for query_name in querys_locs.keys():
            qsite_motif = _get_5mer(query_name, querys_locs[query_name][tlocs_inter[i]], contigname2seq)
            if qsite_motif is None or qsite_motif != tsite_motif:
                is_keep = False
                break
        if is_keep:
            tlocs_inter_tmp.append(tlocs_inter[i])
    tlocs_inter = tlocs_inter_tmp

    if len(tlocs_inter) > 0:
        wstr = ""
        wstr += "{}\t{}\n".format("_".join([target_name, target_strand]),
                                  "\t".join(list(map(str, tlocs_inter))))
        for query_name in querys_locs.keys():
            qlocs = []
            for tloc in tlocs_inter:
                qlocs.append(querys_locs[query_name][tloc])
            wstr += "{}\t{}\n".format(query_name, "\t".join(list(map(str, qlocs))))
        return wstr
    else:
        return ""


# repeat analysis ====================================================
def read_regions(region_file):
    regions = set()
    with open(region_file, "r") as rf:
        for line in rf:
            if line.startswith("##"):
                continue
            words = line.strip().split("\t")
            chrom, start, end = words[0], int(words[1]), int(words[2])
            regions.add((chrom, start, end))
    return regions


def alignment_of_a_group(regions, contigs):
    # 1-based, [a, b]
    regions_len = []
    for region in regions:
        chrom, start, end = region
        regions_len.append((end - start + 1, chrom, start, end))
    regions_len = sorted(regions_len)
    region_ref = regions_len[-1]
    region_querys = regions_len[:-1]

    # temp file
    uuid_tmp = str(uuid.uuid1())
    ref_fa_fp1 = "/tmp" + "/" + "_".join(list(map(str, list(region_ref[1:])))) + "." + uuid_tmp + ".ref+.fa"
    ref_fa_fp2 = "/tmp" + "/" + "_".join(list(map(str, list(region_ref[1:])))) + "." + uuid_tmp + ".ref-.fa"
    with open(ref_fa_fp1, "w") as wf:
        seqlen, chrom, start, end = region_ref
        seq_ref = contigs[chrom][(start-1):end]
        wf.write(">{}\n".format("_".join([chrom, str(start), str(end)])))
        wf.write(seq_ref + "\n")
    with open(ref_fa_fp2, "w") as wf:
        seqlen, chrom, start, end = region_ref
        seq_ref = complement_seq(contigs[chrom][(start-1):end])
        wf.write(">{}\n".format("_".join([chrom, str(start), str(end)])))
        wf.write(seq_ref + "\n")
    query_fa_fp = "/tmp" + "/" + "_".join(list(map(str, list(region_ref[1:])))) + "." + uuid_tmp + ".query.fa"
    with open(query_fa_fp, "w") as wf:
        for region_q in region_querys:
            seqlen, chrom, start, end = region_q
            seq_q = contigs[chrom][(start-1):end]
            wf.write(">{}\n".format("_".join([chrom, str(start), str(end)])))
            wf.write(seq_q + "\n")

    paffp_fwd = _align_minimap2_asm(ref_fa_fp1, query_fa_fp)
    paffp_com = _align_minimap2_asm(ref_fa_fp2, query_fa_fp)
    os.remove(ref_fa_fp1)
    os.remove(ref_fa_fp2)
    os.remove(query_fa_fp)
    return paffp_fwd, paffp_com


def get_match_locs_of_a_group(repeat_file, contigs, motifs, mod_loc):
    repeat_regions = read_regions(repeat_file)
    paf_fwd, paf_com = alignment_of_a_group(repeat_regions, contigs)
    print("==g1")
    wstr_fwd = parse_paf(paf_fwd, contigs, "+", motifs, mod_loc)
    print("==g2")
    wstr_com = parse_paf(paf_com, contigs, "-", motifs, mod_loc)
    os.remove(paf_fwd)
    os.remove(paf_com)
    return wstr_fwd, wstr_com


def get_sites_of_seled_region(args):

    dnaref = DNAReference(args.ref_fa)
    contig2seq = dnaref.getcontigs()
    del dnaref

    if not os.path.exists(args.mum_region):
        raise ValueError("--mum_region not right")
    outdir = args.mum_region.rstrip("/") + ".sites." + args.motifs
    # outdir = args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        for file in os.listdir(outdir):
            os.remove(outdir + "/" + file)
    f_cnt = 0
    valid_cnt1, valid_cnt2 = 0, 0
    for repeat_file in os.listdir(args.mum_region):
        print("=={}: {}".format(f_cnt, repeat_file))
        fname, fext = os.path.splitext(repeat_file)
        repeat_path = "/".join([args.mum_region, repeat_file])
        wstr_fwd, wstr_com = get_match_locs_of_a_group(repeat_path, contig2seq, args.motifs, args.mod_loc)
        valid_cnt = 0
        if wstr_fwd != "":
            valid_cnt += 1
            outfp = outdir + "/" + fname + ".g1.txt"
            with open(outfp, "w") as wf:
                wf.write(wstr_fwd)
        if wstr_com != "":
            valid_cnt += 1
            outfp = outdir + "/" + fname + ".g2.txt"
            with open(outfp, "w") as wf:
                wf.write(wstr_com)
        if valid_cnt == 2:
            valid_cnt2 += 1
        elif valid_cnt == 1:
            valid_cnt1 += 1
        f_cnt += 1
    print("\nrepeat_file_num: {}; has_motif_in_both_strands: {}, "
          "has_motif_in_single_strands: {}".format(f_cnt, valid_cnt2, valid_cnt1))


def main():
    parser = argparse.ArgumentParser("input: mummer region file (chrom, start(1-based), end), genome reference")
    parser.add_argument("--mum_region", type=str, required=True, help="dir")
    parser.add_argument("--ref_fa", type=str, required=True, help="")
    parser.add_argument("--motifs", action="store", type=str,
                        required=False, default='CG',
                        help='motif seq to be extracted, default: CG. '
                             'can be multi motifs splited by comma '
                             '(no space allowed in the input str), '
                             'or use IUPAC alphabet, '
                             'the mod_loc of all motifs must be '
                             'the same')
    parser.add_argument("--mod_loc", action="store", type=int, required=False, default=0,
                        help='0-based location of the targeted base in the motif, default 0')
    # parser.add_argument("--outdir", type=str, required=True, help="")

    args = parser.parse_args()
    get_sites_of_seled_region(args)


if __name__ == '__main__':
    main()
