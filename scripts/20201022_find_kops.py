#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 10:05:49 2019

@author: tomoki
"""

import re
import codecs
import csv, operator
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation

with open("/Users/tomoki/github/kops/analysis/20201022_target_kops.csv", "w") as csv_file:
    fieldnames = ["sequence", "start", "end", "direction"]
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    writer.writeheader()

for record in SeqIO.parse("/Users/tomoki/github/kops/data/GCF_000010245.2_ASM1024v1_genomic.gbff", "genbank"):
    seq = record.seq # Add sequence of E.coli genome

    nlist = ["A", "T", "G", "C"]
    ori = 3677235
    ter = 1561332
    head = 0
    tail = 4646332
    
    for n in nlist:
        dna = "GGG" + n + "AGGG"
        kops = Seq(dna, generic_dna)
        re_kops = kops[::-1]
        comp_kops = re_kops.complement()

        iterator = re.finditer(str(kops), str(seq))
        for match in iterator:
            target = match.group()
            start = match.start()
            end = match.end()
            if ori <= int(match.start()) < tail or head <= int(match.start()) < ter:
                output = str(target) + "," +str(match.start()) + "," + str(match.end()) + "," + str("leading") + "," + seq[start-75:start+75]
                print (output)
            else:
                output = str(target) + "," +str(match.start()) + "," + str(match.end()) + "," + str("lagging") + "," + seq[start-75:start+75]
                print (output)
            with open("/Users/tomoki/github/kops/analysis/20201022_target_kops.csv", "a") as f:
                print (output, file=f)
            
        iterator = re.finditer(str(comp_kops), str(seq))
        for match in iterator:
            target = match.group()
            start = match.start()
            end = match.end()
            if ter <= int(match.start()) < ori:
                output = str(target) + "," +str(match.start()) + "," + str(match.end()) + "," + str("lagging") + "," + seq[start-75:start+75]
                print (output)
            else:
                output = str(target) + "," +str(match.start()) + "," + str(match.end()) + "," + str("leading") + "," + seq[start-75:start+75]
                print (output)
            with open("/Users/tomoki/github/kops/analysis/20201022_target_kops.csv", "a") as f:
                print (output, file=f)
        
    else:
        pass
    
"""
Reference:
    
    https://qiita.com/wanwanland/items/ce272419dde2f95cdabc
    https://python.keicode.com/lang/regular-expression-finditer.php
    
"""