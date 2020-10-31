#! /usr/bin/python
import argparse

sep = "\t"


def extract_posinfo(feafile):
    poses = set()
    with open(feafile, "r") as rf:
        for line in rf:
            words = line.strip().split(sep)
            poses.add(sep.join([words[0], words[1], words[2]]))
    return poses


def write_poses(poses, wfile):
    poses = sorted(list(poses))
    with open(wfile, "w") as wf:
        for pos in poses:
            wf.write(pos + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--feafile", type=str, required=True, help="")

    args = parser.parse_args()

    wfile = args.feafile + ".poses"
    poses = extract_posinfo(args.feafile)
    write_poses(poses, wfile)


if __name__ == '__main__':
    main()
