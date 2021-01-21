#!/usr/bin/env bash

# reads length
cat fqfile | awk '{{if(NR%4==2) print length($1)}}' > fqlenfile

awk '{count++; bases += $1} END{print bases/count}' fqlenfile

# python
import sys
read_lens = []
with open(sys.argv[1]) as rf:
    for line in rf:
        read_lens.append(int(line.strip()))
print("no.: {}, avg len: {}".format(len(read_lens), sum(read_lens)/float(len(read_lens))))

# extract fq
poretools fastq --type fwd --group 001 fast5_pass2.single/10x.1 > fast5_pass2.single.10x.1.fastq &