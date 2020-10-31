#! /usr/bin/python
import argparse


def cal_genome_methy_level(bs_rmet_file):
    met_total, cov_total = 0, 0
    with open(bs_rmet_file, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            met_total += int(words[4])
            cov_total += int(words[6])
    return met_total / float(cov_total)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bs_rmet", type=str, action="append", required=True, help="bs rmet format")

    args = parser.parse_args()

    for bs_rmet in args.bs_rmet:
        print("{}: {}".format(bs_rmet, round(100 * cal_genome_methy_level(bs_rmet), 1)))


if __name__ == '__main__':
    main()
