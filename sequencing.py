"""
pcDNA3のCMVプロモーターとbghpolyAに挟まれた領域のシーケンシング結果を抽出、
リストと照らし合わせて一致するかを判断するプログラム。
list.txtを準備し、aliasと組み合わせてseqのみで実行可能なように設定した。
"""

import sys
import glob
import re
file = sys.argv

#第一変数にリスト名
#リストは配列のみの一次元、あるいはサンプル名と配列のタブ区切り二次元形式。

def extractionA(a): #いつものやつにファイル形式の判定を行う機能を追加した。タブの有無で一か二次元を判定。modeはSingleかW。
	with open(a) as fi:
		b = fi.readlines()
		if "\t" in b[0]:
			c = [i.rstrip().split("\t") for i in b]
			mode = "W"
		else:
			c = [i.rstrip() for i in b]
			mode = "S"
		return c, mode

def cutseq(a,b,c): #.seqファイルからbとcで指定された配列の間を切り出す。
	with open(a) as fi:
		d = fi.read().replace("\n","")
		e = re.sub(".*" + b,"",re.sub(c + ".*", "", d))
		return e

#切り出す両端の配列を指定。
Fw = "AAGCTTGCTAGCCACC"
Rev = "GGATCCGGCGGCTCCGGAGGAGGTACC"


ref, mode = extractionA(file[1])
seqdata = glob.glob("*.seq")
seqdata.sort()

atari, hazure = 0,0

for i in range(len(ref)):
	reference = ref[i]
	if mode == "S":
		refname = ""
		reference = ref[i]
	elif mode == "W":
		refname = ref[i][0]
		reference = ref[i][1]

	result = cutseq(seqdata[i],Fw,Rev)

	print(str(i+1) + "\t", end="")
	if reference == result:
		print("あたり", end = "")
		atari+=1
	else:
		print("ハズレ", end = "")
		hazure+=1
	print("\t" + refname)
	print("Reference\t" + reference)
	print("Seq_Result\t" + result)
	print()
print("correct seq\t" + str(atari))
print("wrong seq\t" + str(hazure))
print("")
if hazure > atari:
	print("お前やる気あんの？")
elif atari >= hazure and atari < hazure * 2:
	print("もうちょっと頑張れよ")
else:
	print("自惚れるのは早いぞ")
print("")
