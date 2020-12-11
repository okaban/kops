import re
import codecs
import sys
import csv, operator
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation

args = sys.argv

for rec in SeqIO.parse(args[1], 'genbank'):
    source = rec.annotations['source']
    name = rec.name
    length = len(rec.seq)

for record in SeqIO.parse(args[1], "genbank"):
    seq = record.seq
 
    import pandas as pd #pandasをインポート
    df = pd.read_csv(args[2], encoding='shift_jis')    #CSVファイルを読み込む
    pre_ori = df.iat[10,1]  #oriのポジジョンを読み込む
    pre_ter = df.iat[10, 3] #terのポジジョンを読み込む

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

        #リーディング鎖のKOPSをカウント
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

        #ラギング鎖のKOPSをカウント
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

#print ("Predicted origin is " + str(ori))
#print ("Predicted terminus is " + str(ter))
#print ("Leading count is " + str(lead_count))
#print ("Lagging count is " + str(lagg_count))

kops_skew = lead_count / (lead_count + lagg_count)
#print ("KOPS skew is " + str(kops_skew))

print (str(name)+"\t"+str(source)+"\t"+str(length)+"\t"+str(ori)+"\t"+str(ter)+"\t"+str(lead_count)+"\t"+str(lagg_count)+"\t"+str(kops_skew))