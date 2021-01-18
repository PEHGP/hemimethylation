#!/usr/bin/env python2.7
#coding=utf-8
#直接过滤掉了小于阈值的点
import sys,gzip
#bismark_methylation_extractor --CX
f=sys.argv[1]#RnaseH1a2a_pe.deduplicated.CX_report.txt.gz
cov=sys.argv[2]#4
prefix=sys.argv[3]
d={}
for x in gzip.open(f):
	x=x.rstrip()
	l=x.split("\t")
	if not l[5] in d:
		d[l[5]]={}
	count_methylated=float(l[3])
	count_nonmethylated=float(l[4])
	if not l[0] in d[l[5]]:
		d[l[5]][l[0]]=[]
	if count_methylated+count_nonmethylated<float(cov):
		#d[l[5]][l[0]].append((int(l[1]),"NaN",l[2]))
		continue
	d[l[5]][l[0]].append((int(l[1]),count_methylated/(count_methylated+count_nonmethylated),l[2]))
for co in d:
	Fr=open(prefix+"_"+co+".wig","w")
	for ch in d[co]:
		Fr.write("variableStep chrom=%s\n"%ch)
		for m in sorted(d[co][ch]):
			if m[1]==0:
				s=str(int(m[1]))
			elif m[2]=="-" and m[1]!="NaN":
				s=str(m[1]*-1)
			else:
				s=str(m[1])
			Fr.write(str(m[0])+"\t"+s+"\n")
	Fr.close()
