import re
import codecs
import csv, operator
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation

with open("/Users/tomoki/github/kops/analysis/kops_position_20201210.csv", "w") as csv_file:
    fieldnames = ["sequence", "start", "end", "direction"]
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    writer.writeheader()

for record in SeqIO.parse("/Users/tomoki/my_projects/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.gbff", "genbank"):
    seq = record.seq
 
    import pandas as pd #pandasをインポート
    df = pd.read_csv('/Users/tomoki/my_projects/kops/skew/scripts/result.csv', encoding='shift_jis')    #CSVファイルを読み込む
    pre_ori = df.iat[10,1]  #oriのポジジョンを読み込む
    pre_ter = df.iat[10, 3] #terのポジジョンを読み込む

    nlist = ["A", "T", "G", "C"]
    ori = int(pre_ori)
    ter = int(pre_ter)
    head = 0
    tail = 4646332
    
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
                #print (output)
                lead_count += 1

        iterator = re.finditer(str(re_comp_kops), str(seq))
        for match in iterator:
            target = match.group()
            start = match.start()
            end = match.end()
            if ter <= int(match.start()) < ori:    #検索範囲をラギング鎖（ter-ori）に指定
                output = str(target) + "," +str(match.start()+1) + "," + str(match.end()) + "," + str("leading") + "," + seq[start-75:start+75]
                #print (output)
                lead_count += 1

        #ラギング鎖のKOPSをカウント
        iterator = re.finditer(str(re_comp_kops), str(seq))
        for match in iterator:
            target = match.group()
            start = match.start()
            end = match.end()
            if ori <= int(match.start()) < tail or head <= int(match.start()) < ter:    #検索範囲をリーディング鎖（ori-tailとhead-ter）に指定
                output = str(target) + "," +str(match.start()+1) + "," + str(match.end()) + "," + str("lagging") + "," + seq[start-75:start+75]
                #print (output)
                lagg_count += 1

        iterator = re.finditer(str(kops), str(seq))
        for match in iterator:
            target = match.group()
            start = match.start()
            end = match.end()
            if ter <= int(match.start()) < ori:    #検索範囲をラギング鎖（ter-ori）に指定
                output = str(target) + "," +str(match.start()+1) + "," + str(match.end()) + "," + str("lagging") + "," + seq[start-75:start+75]
                #print (output)
                lagg_count += 1

    else:
        pass

print ("Predicted origin is " + str(ori))
print ("Predicted terminus is " + str(ter))
print ("Leading count is " + str(lead_count))
print ("Lagging count is " + str(lagg_count))

kops_skew = lead_count / (lead_count + lagg_count)
print ("KOPS skew is " + str(kops_skew))