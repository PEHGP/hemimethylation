#!/usr/bin/env python2.7
#coding=utf-8
#直接过滤掉了小于阈值的点
import sys,gzip,os
def GetCGWig(f,cov,prefix):
	Fr=open(prefix+"_CG.wig","w")
	ch=[]
	lines=os.popen("zgrep -P '\tCG\t' %s"%f).readlines()
	for i,x in enumerate(lines):
		x=x.rstrip()
		l=x.split("\t")
		if not l[0] in ch:
			ch.append(l[0])
			Fr.write("variableStep chrom=%s\n"%l[0])
		if l[2]=="-":
			continue
		rev=lines[i+1].rstrip().split("\t")
		if rev[2]!="-" or int(rev[1])!=int(l[1])+1:
			print("error")
			print(l)
			print(rev)
			sys.exit()
		fwd_m=float(l[3])
		fwd_um=float(l[4])
		rev_m=float(rev[3])
		rev_um=float(rev[4])
		if fwd_m+fwd_um<=float(cov) and rev_m+rev_um<=float(cov):
			continue
		if fwd_m+fwd_um==0:
			Fr.write(l[1]+"\t"+str(0)+"\n")
		else:
			Fr.write(l[1]+"\t"+str(fwd_m/(fwd_m+fwd_um))+"\n")
		if rev_m+rev_um==0:
			Fr.write(rev[1]+"\t"+str(0)+"\n")
		else:
			Fr.write(rev[1]+"\t"+str(rev_m/(rev_m+rev_um)*-1)+"\n")
	Fr.close()
def GetCHGWig(f,cov,prefix):
	Fr=open(prefix+"_CHG.wig","w")
	ch=[]
	lines=os.popen("zgrep -P '\tCHG\t' %s"%f).readlines()
	for i,x in enumerate(lines):
		x=x.rstrip()
		l=x.split("\t")
		if not l[0] in ch:
			ch.append(l[0])
			Fr.write("variableStep chrom=%s\n"%l[0])
		if l[2]=="-":
			continue
		if l[-1]=="CCG":
			if float(l[3])+float(l[4])<=float(cov):
				continue
			else:
				Fr.write(l[1]+"\t"+str(float(l[3])/(float(l[3])+float(l[4])))+"\n")
				continue
		rev=lines[i+1].rstrip().split("\t")
		if rev[2]!="-" or int(rev[1])!=int(l[1])+2:
			print("error")
			print(l)
			print(rev)
			sys.exit()
		fwd_m=float(l[3])
		fwd_um=float(l[4])
		rev_m=float(rev[3])
		rev_um=float(rev[4])
		if fwd_m+fwd_um<=float(cov) and rev_m+rev_um<=float(cov):
			continue
		if fwd_m+fwd_um==0:
			Fr.write(l[1]+"\t"+str(0)+"\n")	
		else:
			Fr.write(l[1]+"\t"+str(fwd_m/(fwd_m+fwd_um))+"\n")
		if rev_m+rev_um==0:
			Fr.write(rev[1]+"\t"+str(0)+"\n")
		else:
			Fr.write(rev[1]+"\t"+str(rev_m/(rev_m+rev_um)*-1)+"\n")
	Fr.close()
def GetCHHWig(f,cov,prefix):
	Fr=open(prefix+"_CHH.wig","w")
	ch=[]
	lines=os.popen("zgrep -P '\tCHH\t' %s"%f).readlines()
	for i,x in enumerate(lines):
		x=x.rstrip()
		l=x.split("\t")
		if not l[0] in ch:
			ch.append(l[0])
			Fr.write("variableStep chrom=%s\n"%l[0])
		if float(l[3])+float(l[4])<=float(cov):
				continue
		Fr.write(l[1]+"\t"+str(float(l[3])/(float(l[3])+float(l[4])))+"\n")
	Fr.close()
if __name__ == '__main__':
	#bismark_methylation_extractor --CX
	f=sys.argv[1]#RnaseH1a2a_pe.deduplicated.CX_report.txt.gz
	cov=sys.argv[2]#4
	prefix=sys.argv[3]
	GetCGWig(f,cov,prefix)
	GetCHGWig(f,cov,prefix)
	GetCHHWig(f,cov,prefix)
