import argparse
import os


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


def cut_genome_N(contignames, contigs):
    contignames_new, contigs_new = [], dict()
    print("==================== cut N")
    for cname in contignames:
        contigseq = contigs[cname]
        start = 0
        end = 0
        for idx in range(1, len(contigseq)):
            lastbase = contigseq[idx-1]
            currbase = contigseq[idx]
            if lastbase == "N" and currbase == "N":
                continue
            elif lastbase == "N" and currbase != "N":
                start = idx
                end = idx
            elif lastbase != "N" and currbase != "N":
                end = idx
            elif lastbase != "N" and currbase == "N":
                # subchr_name = "_".join([cname, str(start), str(end+1)])
                subchr_name = (cname, start, end+1)
                subchr_seq = contigseq[start:(end+1)]
                contignames_new.append(subchr_name)
                contigs_new[subchr_name] = subchr_seq
                print(subchr_name)
        if end == len(contigseq) - 1:
            # subchr_name = "_".join([cname, str(start), str(end + 1)])
            subchr_name = (cname, start, end + 1)
            subchr_seq = contigseq[start:(end + 1)]
            contignames_new.append(subchr_name)
            contigs_new[subchr_name] = subchr_seq
            print(subchr_name)
    # for cname in contignames_new:
    #     print(cname)
    return contignames_new, contigs_new


def combine_subchrseqs(contignames, contigs):
    # contigsnames must be sorted
    contignames_new, contigs_new = [], dict()
    last_cname, last_s, last_e = contignames[0]
    last_seq = contigs[(last_cname, last_s, last_e)]
    print("============== recombine contig")
    for idx in range(1, len(contignames)):
        curr_cname, curr_s, curr_e = contignames[idx]
        if curr_cname != last_cname or curr_s - last_e >= 20:
            contignames_new.append((last_cname, last_s, last_e))
            contigs_new[(last_cname, last_s, last_e)] = last_seq
            print(last_cname, last_s, last_e)

            last_cname, last_s, last_e = curr_cname, curr_s, curr_e
            last_seq = contigs[(last_cname, last_s, last_e)]
        else:
            last_seq = last_seq + "N" * (curr_s - last_e) + contigs[contignames[idx]]
            last_e = curr_e
    contignames_new.append((last_cname, last_s, last_e))
    contigs_new[(last_cname, last_s, last_e)] = last_seq
    print(last_cname, last_s, last_e)
    return contignames_new, contigs_new


def write_contig_to_fa(contignames, contigs, wfile):
    wf = open(wfile, "w")
    for cname in contignames:
        wf.write(">{}\n".format("_".join(list(map(str, cname)))))
        wf.write(contigs[cname] + "\n")
    wf.close()


def main():
    parser = argparse.ArgumentParser("remove Ns of genome")
    parser.add_argument("--genome", type=str, required=True, help="")

    args = parser.parse_args()
    dnaref = DNAReference(args.genome)
    contignames_new, contigs_new = cut_genome_N(dnaref.getcontignames(), dnaref.getcontigs())
    contignames_new, contigs_new = combine_subchrseqs(contignames_new, contigs_new)
    fname, fext = os.path.splitext(args.genome)
    wfile = fname + ".rmN_cut" + fext
    write_contig_to_fa(contignames_new, contigs_new, wfile)


if __name__ == '__main__':
    main()
