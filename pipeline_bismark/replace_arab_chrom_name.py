import sys
import os

chromname_map_arab = {"NC_003070.9": "Chr1",
                      "NC_003071.7": "Chr2",
                      "NC_003074.8": "Chr3",
                      "NC_003075.7": "Chr4",
                      "NC_003076.8": "Chr5",
                      "NC_037304.1": "ChrM",
                      "NC_000932.1": "ChrC"}


def main(argv):
    rfile = argv[1]
    fname, fext = os.path.splitext(rfile)
    wfile = fname + ".chgname" + fext
    wf = open(wfile, "w")
    with open(rfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            try:
                words[0] = chromname_map_arab[words[0]]
                wf.write("\t".join(words) + "\n")
            except KeyError:
                continue
    wf.close()


if __name__ == '__main__':
    main(sys.argv)
