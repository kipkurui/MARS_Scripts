#!/usr/bin/python
# Philip Machanick

from __future__ import print_function
import sys
from math import log
import os

found = 0
row = 0
n_rows = 0
tn_rows = 0
entropy = 0
total_entropy = 0
motifs = 0
name = ""
if len(sys.argv) < 3:
    print("Usage: python motif-IC.py input.meme")
    sys.exit(0)
n = 0
motif_file = sys.argv[1]
test = motif_file
motif_file2 =test
rwa_file = sys.argv[2]
tf_name = sys.argv[3]
raw_dict = {}
with open(rwa_file) as raw_in:
    for line in raw_in:
        raw_dict[line.split()[0]] = line.split()[-1]
out = "Motif_name\tMotif_IC\tMotif_IC\tMotif_length\tMotif_score\tMotif_logo"
print(out)
with open(motif_file, "r") as motif_file:
    for line in motif_file:
        words = line.split()
        if found == 0:
            if line.startswith("MOTIF"):
                # allow for motifs without an alternative name
                if len(words) < 3:
                    words.append("")
                name = (words[1])
                found = 1
                motifs += motifs
                entropy = 0
                continue
        if found == 1:
            if line.startswith("letter-probability"):
                n_rows = int((line.split("w="))[1].split()[0])
                found = 2
            continue
        if found == 2:
            if line == "\n":
                continue
            else:
                check = 0
            for val in words:
                if float(val) > 0:
                    check += float(val) * log(float(val))/log(2.0)
                    entropy += float(val) * log(float(val))/log(2.0)
            row += 1
            if row >= n_rows:
                v = 2*n_rows+entropy
                out = "%s\t%f\t%f\t%i\t%f\t<img src='/static/files/%s/motifs/%s.png' alt='My image' class='img-responsive'/>"\
                      % (name, v, (v/n_rows), n_rows, float(raw_dict[name]), tf_name, name)
                print(out)
                #n+= 1
                #print(n)
                found = 0
                row = 0
                total_entropy += (v/n_rows)

mot_list = []
os.system("mkdir -p /home/kipkurui/Project/MARS/MATOM/static/files/%s/motifs" % tf_name)
with open(test) as meme_in:
    for line in meme_in:
        if line.startswith("MOTIF"):
            mot_list.append(line.split()[1])
fold = "/home/kipkurui/Project/MARS/MATOM/static/files/%s/motifs" % tf_name
for i in mot_list:
    os.system("ceqlogo -i%s %s -f PNG -h3 -w4 -o %s/%s.png" % (i, motif_file2, fold, i))