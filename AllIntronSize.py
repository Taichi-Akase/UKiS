#送信の都合上拡張子を.txtで送っているので.pyにリネームしてください
#元となるデータ(RefSeq.csv)はNCBIのTable Browserより獲得し、事前に同名の遺伝子(txStart/txEndを参考に転写領域が重複しているもの)を取り除いています
#Refseq.csvが保存された場所のパスを下記os.chdirに記してもらえれば準備完了です
#https://note.nkmk.me/python-os-getcwd-chdir/
#あとはpandasとmatplotlibが導入されたpython3にて実行が可能です

#AllIntronSize.py:56: FutureWarning: The frame.append method is deprecated and will be removed from pandas in a future version. Use pandas.concat instead.
#上記警告文のため59行のIntron_All.appendをpd.concatに変更。20220522_TA

import pandas as pd
import os
import math
import numpy as np
import matplotlib.pyplot as plt # グラフ描画
#os.chdir(r"") #RefSeqデータ格納場所指定
os.chdir(os.getcwd())

#RefSeqデータ由来の各遺伝子のエキソン位置から全イントロン長を算出
Intron=pd.read_csv("RefSeq.csv")
#RefseqのExon位置データからIntron位置データを出す
Intron['IntronStarts'] = Intron['exonEnds'].str.split(',')
Intron['IntronEnds'] = Intron['exonStarts'].str.split(',')
Intron['IntronSize'] = ""
for i in range(len(Intron)):
    Starts=""
    Ends=""
    for j in range(Intron['exonCount'][i]-1):
        Starts=Starts+str(int(Intron['IntronStarts'][i][j])+1)+","
        Ends=Ends+str(int(Intron['IntronEnds'][i][j+1])-1)+","
    Intron['IntronStarts'][i]=Starts
    Intron['IntronEnds'][i]=Ends
Intron['StartList']=Intron['IntronStarts'].str.split(',')
Intron['EndList']=Intron['IntronEnds'].str.split(',')

#イントロン長を導出する
for i in range(len(Intron)):
    allSize=''
    for j in range(Intron['exonCount'][i]-1):
        Size=int(Intron['EndList'][i][j])-int(Intron['StartList'][i][j])+1
        if j==0:
            allSize=str(Size)
        else:
            allSize=allSize+','+str(Size)
        Size=0
    Intron['IntronSize'][i]=allSize+','

#遺伝子ストランドから第一イントロン長のみ抽出
Intron_First=Intron
Intron_First['lengthList'] = Intron_First['IntronSize'].str.split(',')
Intron_First['firstIntronSize'] = [x[0] if Intron_First['strand'][count]=='+' else x[-2] for (count,x) in enumerate(Intron_First['lengthList'])]
Intron_First=Intron_First[Intron_First['firstIntronSize']!=""]

#全イントロン長をリストアップ
data_all= Intron['IntronSize'].str.split(',')
Intron_All = pd.DataFrame(index=[],columns=['IntronSize'])
for i in range(len(data_all)):
    for j in range(len(data_all[i])-1):
        series = pd.Series([data_all[i][j]], index=Intron_All.columns)
        Intron_All = Intron_All.append(series, ignore_index = True)
        #Intron_All = pd.concat(series, ignore_index = True)
Intron_All=Intron_All[Intron_All['IntronSize']!=""]

Intron_First.to_csv("Intron_first.csv")
Intron_All.to_csv("Intron_all.csv")


# 箱ひげ図作製
data_First=pd.read_csv("Intron_First.csv")
data_All=pd.read_csv("Intron_All.csv")
data=data_First['firstIntronSize']
data_all= data_All['IntronSize']

fig, ax = plt.subplots(figsize=(8, 16))
plt.rcParams["font.size"] = 20
plt.yscale("log") # y軸対数グラフ

bp = ax.boxplot([data,data_all],
                        whis="range", widths=.75,
                        showmeans=True,
                        meanprops=dict(marker='+', markersize=20,markeredgecolor='y')) #外れ値も含む,平均値も記す
TP53=10753 #TP53第一イントロン長
bp = ax.hlines(TP53, 0, 3, "blue", linestyles='dashed')
ax.set_xticklabels(['firstIntron','allIntron'])
# タイトル名
plt.title('IntronSize')
# ラベル名
#plt.xlabel('subject')
plt.ylabel('length [bp]')
# Y軸のメモリの長さ
plt.ylim([-1,3000000])
#plt.grid()
plt.savefig("IntronSize.png")
# 描画
plt.show()
