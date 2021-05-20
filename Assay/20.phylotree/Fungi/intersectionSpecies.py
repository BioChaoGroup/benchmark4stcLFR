import sys
import os
import re
from collections import defaultdict


mydic = {}
test_file = r'../../../Result/phylotree/nature05110/01.data/genus.both.inf.txt'
with open(test_file, 'r') as fh:
    for line in fh.readlines():
        line = line.strip()
        tmp = line.split('\t')
        #spe = tmp[1].split("species__")
        if len(tmp) == 2:
            mydic[tmp[1]] = tmp[0]


iddic = {}
ifa = open('../../../Result/phylotree/nature05110/01.data/species.both.mafft.trimal.fa.reduced', "r")
while True:
    line1 = ifa.readline()
    if line1:
        line1 = line1.split(" ")
        iddic[line1[0]] = line1[1]
    else:
        break
ifa.close()


file = r'../../../Result/phylotree/nature05110/01.data/001_41586_2006_BFnature05110_MOESM1_ESM.txt'
with open(file, 'r') as fa:
    fa.readline()
    for li in fa.readlines():
        li = li.strip()
        li = li.split("\t")
        if li[2] in mydic.keys():
            if mydic[li[2]] in iddic.keys():
                print(">",li[2])
                print(iddic[mydic[li[2]]], end="")