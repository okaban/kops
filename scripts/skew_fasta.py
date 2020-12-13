#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 13 2:36:20 2020

@author: tomoki
"""

#   モジュールをインポート
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import pandas as pd
import csv, operator
import sys
import re

args = sys.argv
df = pd.read_csv(args[2], encoding='shift_jis')    #CSVファイルを読み込む
pre_ori = df.iat[10,1]  #oriのポジジョンを読み込む
pre_ter = df.iat[10, 3] #terのポジジョンを読み込む

#   FASTAファイルを読み込む
records_dict = SeqIO.to_dict(SeqIO.parse(args[1],"fasta"))
for record in SeqIO.parse(args[1], "fasta"):
    seq = record.seq
    description = record.description
    id = record.id
    length = len(record.seq)

    nlist = ["A", "T", "G", "C"]
    ori = int(pre_ori)
    ter = int(pre_ter)
    head = 0
    tail = int(length)

    lead_count = 0
    lagg_count = 0

    for n in nlist:
        dna = "GGG" + n + "AGGG"
        kops = Seq(dna, generic_dna)
        re_kops = kops[::-1]
        comp_kops = kops.complement()
        re_comp_kops = re_kops.complement()

#   リーディング鎖のKOPSをカウント
        iterator = re.finditer(str(kops), str(seq))
        for match in iterator:
            target = match.group()
            start = match.start()
            end = match.end()
            if ori <= int(match.start()) < tail or head <= int(match.start()) < ter:    #検索範囲をリーディング鎖（ori-tailとhead-ter）に指定
                output = str(target) + "," +str(match.start()+1) + "," + str(match.end()) + "," + str("leading") + "," + seq[start-75:start+75]
                lead_count += 1

        iterator = re.finditer(str(re_comp_kops), str(seq))
        for match in iterator:
            target = match.group()
            start = match.start()
            end = match.end()
            if ter <= int(match.start()) < ori:    #検索範囲をラギング鎖（ter-ori）に指定
                output = str(target) + "," +str(match.start()+1) + "," + str(match.end()) + "," + str("leading") + "," + seq[start-75:start+75]
                lead_count += 1

#   ラギング鎖のKOPSをカウント
        iterator = re.finditer(str(re_comp_kops), str(seq))
        for match in iterator:
            target = match.group()
            start = match.start()
            end = match.end()
            if ori <= int(match.start()) < tail or head <= int(match.start()) < ter:    #検索範囲をリーディング鎖（ori-tailとhead-ter）に指定
                output = str(target) + "," +str(match.start()+1) + "," + str(match.end()) + "," + str("lagging") + "," + seq[start-75:start+75]
                lagg_count += 1

        iterator = re.finditer(str(kops), str(seq))
        for match in iterator:
            target = match.group()
            start = match.start()
            end = match.end()
            if ter <= int(match.start()) < ori:    #検索範囲をラギング鎖（ter-ori）に指定
                output = str(target) + "," +str(match.start()+1) + "," + str(match.end()) + "," + str("lagging") + "," + seq[start-75:start+75]
                lagg_count += 1

    else:
        pass

#   KOPS-skewを計算
kops_skew = lead_count / (lead_count + lagg_count)

#   結果をタブ区切りで出力
print (str(id)+"\t"+str(description.lstrip(id))+"\t"+str(length)+"\t"+str(ori)+"\t"+str(ter)+"\t"+str(lead_count)+"\t"+str(lagg_count)+"\t"+str(kops_skew))