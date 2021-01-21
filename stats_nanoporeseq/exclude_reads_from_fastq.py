#!/usr/bin/python
import argparse


def _exclude_read_ids_from_fastq(fqfile, wfqfile, readidfile):
    readids = set()
    with open(readidfile, "r") as rf:
        for line in rf:
            readids.add(line.strip())
    wf = open(wfqfile, "w")
    with open(fqfile, 'r') as rf:
        lidx = 0
        iswrite = True
        for line in rf:
            if lidx % 4 == 0:
                assert(line.startswith("@"))
                readid = line.strip().split(" ")[0][1:]
                if readid in readids:
                    iswrite = False
                else:
                    iswrite = True
                    wf.write(line)
            else:
                if iswrite:
                    wf.write(line)
            lidx += 1
    wf.flush()
    wf.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fqfile", type=str, help="")
    parser.add_argument("--ridfile", type=str, help="readids to be excluded in fqfile")
    parser.add_argument("--wfqfile", type=str, default=None,
                        help="output fqfile")

    args = parser.parse_args()

    _exclude_read_ids_from_fastq(args.fqfile, args.wfqfile, args.ridfile)


if __name__ == '__main__':
    main()
