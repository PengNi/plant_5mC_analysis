#!/usr/bin/python
import argparse
import os

same_region_cf = 100
# identity_cf = 0.95
over_region_cf = 100
# over_region_ratio = 0.3
# len_cf = 9000


class MumRepeat:
    # [start, end], 0-based
    def __init__(self, fields):
        self._refstart = int(fields[0])
        self._refend = int(fields[1])
        self._querystart = int(fields[2])
        self._queryend = int(fields[3])
        self._ref_regionlen = int(fields[4])
        self._query_regionlen = int(fields[5])
        self._identity = float(fields[6])
        self._refid = fields[11]
        self._queryid = fields[12]

        self._ref_strand = "+"
        self._query_strand = "+"
        if self._querystart > self._queryend:
            self._query_strand = "-"
            qstart = self._querystart
            self._querystart = self._queryend
            self._queryend = qstart

    def get_repeats_ref2ref_nosplit(self, len_cf, identity_cf, contigset, contig_prefix):
        if contigset is not None:
            if (self._refid not in contigset) or (self._queryid not in contigset):
                return []
        if contig_prefix is not None:
            if (not str(self._refid).startswith(contig_prefix)) or (not str(self._queryid).startswith(contig_prefix)):
                return []

        if self._identity < identity_cf or (self._ref_regionlen < len_cf and self._query_regionlen < len_cf):
            return []
        if self._refid == self._queryid:
            # if abs(self._refstart - self._querystart) < same_region_cf * (1-over_region_ratio) and \
            #                 abs(self._refend - self._queryend) < same_region_cf * (1-over_region_ratio):
            #     return []
            # minend = min(self._refend, self._queryend)
            # maxstart = max(self._refstart, self._querystart)
            # over_region_cf = (self._refend - self._refstart) * over_region_ratio
            # if maxstart <= minend - over_region_cf:
            #     print(self._refid, self._refstart, self._refend, self._querystart, self._queryend)
            #     return [(self._refid, min(self._refstart, self._querystart),
            #             max(self._refend, self._queryend))]
            if self._refstart == self._querystart and self._refend == self._queryend:
                return []
        if (self._refid, self._refstart, self._refend) < (self._queryid, self._querystart, self._queryend):
            return [(self._refid, self._refstart, self._refend, self._queryid, self._querystart, self._queryend)]
        else:
            return [(self._queryid, self._querystart, self._queryend, self._refid, self._refstart, self._refend)]

    def get_repeats_ref2ref(self, len_cf, identity_cf):
        """
        find interspected repeats/tandem repeats, two repeat sequences no overlap
        :param len_cf:
        :param identity_cf:
        :return:
        """
        if self._identity < identity_cf or self._ref_regionlen < len_cf:
            return []
        if self._refid == self._queryid:
            # if abs(self._refstart - self._querystart) < same_region_cf and \
            #                 abs(self._refend - self._queryend) < same_region_cf:
            #     return []
            minend = min(self._refend, self._queryend)
            maxstart = max(self._refstart, self._querystart)
            # if maxstart <= minend - over_region_cf:
            if maxstart < minend:
                print(self._refid, self._refstart, self._refend, self._querystart, self._queryend)
                return [(self._refid, min(self._refstart, self._querystart),
                        max(self._refend, self._queryend))]
        return [(self._refid, self._refstart, self._refend),
                (self._queryid, self._querystart, self._queryend)]

    def get_query_repeats_with_ref(self, len_cf, identity_cf):
        if self._identity < identity_cf or self._ref_regionlen < len_cf:
            return None
        query_region = (self._queryid, self._querystart, self._queryend)
        ref_region = (self._refid, self._refstart, self._refend)
        return query_region, ref_region


def _print_repeats(repeats, wfile):
    wf = open(wfile, "w")
    for repeat in repeats:
        wf.write("\t".join(list(map(str, repeat))) + "\n")
    wf.close()


def extract_repeats_mummer(args):
    ref_mumrepeats = set()

    contigset = set(args.contig_names.strip().split(",")) if args.contig_names is not None else None
    contig_prefix = args.contig_prefix
    with open(args.crds_ref_file, 'r') as rf:
        for line in rf:
            words = line.strip().split("\t")
            if len(words) < 13:
                continue
            mumrepeat = MumRepeat(words)
            for mmr in mumrepeat.get_repeats_ref2ref_nosplit(args.len_cf, args.identity_cf, contigset, contig_prefix):
                ref_mumrepeats.add(mmr)
    ref_mumrepeats = sorted(list(ref_mumrepeats))
    # ref_mumrepeats = _unique_repeats(ref_mumrepeats)
    _print_repeats(ref_mumrepeats, args.wfile)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--crds_ref_file", type=str, required=True, help="nucmer show-coords ref2ref file, "
                                                                         "or contig2contig file")
    parser.add_argument("--len_cf", type=int, required=False, default=5000,
                        help="repeat len cut off")
    parser.add_argument("--identity_cf", type=float, required=False, default=95,
                        help="identity cut off, 0-100, default 95")
    parser.add_argument("--wfile", type=str, required=True,
                        help="")
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12

    args = parser.parse_args()
    extract_repeats_mummer(args)


if __name__ == '__main__':
    main()
