#!/usr/bin/python
import argparse
import os
import fnmatch


def _get_read_ids_from_fastq(fqfile):
    rids = set()
    with open(fqfile, 'r') as rf:
        lidx = 0
        for line in rf:
            if lidx % 4 == 0:
                assert(line.startswith("@"))
                rids.add(line.strip().split(" ")[0][1:])
            lidx += 1
    return rids


def move_f5s_2_newdir(args):
    readids = _get_read_ids_from_fastq(args.fqfile)
    print("{} reads selected..".format(len(readids)))

    f5paths = []
    for root, dirnames, filenames in os.walk(args.f5dir):
        for filename in fnmatch.filter(filenames, '*.fast5'):
            fast5_path = os.path.join(root, filename)
            f5paths.append(fast5_path)
    print("{} fast5 files in fast5 dir".format(len(f5paths)))

    f5paths_left = []
    for f5path in f5paths:
        readidtmp = os.path.basename(f5path)[:36]
        if readidtmp in readids:
            f5paths_left.append(f5path)
    print("{} selected reads get f5path mappings".format(len(f5paths_left)))

    ends = []
    for end in range(0, len(f5paths_left), 4000):
        ends.append(end)
    if ends[-1] != len(f5paths_left):
        ends.append(len(f5paths_left))

    for i in range(0, len(ends)-1):
        dstart, dend = ends[i], ends[i+1]

        curr_dir = "/".join([args.opdir, str(i)])
        os.makedirs(curr_dir)

        for f5path in f5paths_left[dstart:dend]:
            os.rename(f5path, curr_dir + "/" + os.path.basename(f5path))

    print("file moving end..")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--f5dir", type=str, help="")
    parser.add_argument("--fqfile", type=str, help="")
    parser.add_argument("--opdir", type=str, default=None,
                        help="output dir")

    args = parser.parse_args()
    move_f5s_2_newdir(args)


if __name__ == '__main__':
    main()
