#! /usr/bin/python
import sys
import os


def main():
    signalfile = sys.argv[1]  # 17.360.file from extract_features.py
    kmerlen = int(sys.argv[2])
    if len(sys.argv) > 3:
        MAX_BASE_NUM = int(sys.argv[3])
    if kmerlen > MAX_BASE_NUM:
        print("kmerlen must <= {}".format(MAX_BASE_NUM))
        return
    fname, fext = os.path.splitext(signalfile)
    wfp = fname + '.' + str(kmerlen) + 'mer.signalnum.txt'

    wf = open(wfp, 'w')
    with open(signalfile, 'r') as rf:
        next(rf)
        count = 0
        for line in rf:
            words = line.strip().split('\t')
            signal_lens = [int(x) for x in words[9].strip().split(',')]
            if kmerlen != MAX_BASE_NUM:
                cut_per_side = (MAX_BASE_NUM - kmerlen) // 2
                signal_lens = signal_lens[cut_per_side:(-cut_per_side)]
            if kmerlen != 1:
                signalnum = sum(signal_lens)
                wf.write(str(signalnum) + '\n')
                count += 1
            else:
                for slen in signal_lens:
                    if slen < 100:
                        wf.write(str(slen) + '\n')
                        count += 1
            if count >= 10000000:
                break
    wf.close()


if __name__ == '__main__':
    main()
