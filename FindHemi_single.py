#!/usr/bin/env python
import sys,pandas,os,collections
import numpy as np
import statsmodels.api as sm
import scipy.stats as stats
def GetCGHemi(f,cutoff):
	r=[]
	lines=os.popen("zgrep -P '\tCG\t' %s"%f).readlines()
	for i,x in enumerate(lines):
		x=x.rstrip()
		l=x.split("\t")
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
		oddsratio, pvalue = stats.fisher_exact([[fwd_m,rev_m],[fwd_um,rev_um]])
		if fwd_m==0:
			fwd_m=epsilon
		fwd_p=fwd_m/(fwd_m+fwd_um)
		if rev_m==0:
			rev_m=epsilon
		rev_p=rev_m/(rev_m+rev_um)
		#fc=np.log2(fwd_p/rev_p)
		diff=abs(fwd_p-rev_p)
		if pvalue<=cutoff:
			r.append((l[0],l[1],pvalue,str(diff)))
		#print(l[0],l[1],[fwd_m,rev_m,fwd_um,rev_um],fwd_p,rev_p,pvalue)
	return r
def GetCHGHemi(f,cutoff):
	r=[]
	lines=os.popen("zgrep -P '\tCHG\t' %s"%f).readlines()
	for i,x in enumerate(lines):
		x=x.rstrip()
		l=x.split("\t")
		if l[2]=="-":
			continue
		if l[-1]=="CCG":
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
		oddsratio, pvalue = stats.fisher_exact([[fwd_m,rev_m],[fwd_um,rev_um]])
		if fwd_m==0:
			fwd_m=epsilon
		fwd_p=fwd_m/(fwd_m+fwd_um)
		if rev_m==0:
			rev_m=epsilon
		rev_p=rev_m/(rev_m+rev_um)
		#fc=np.log2(fwd_p/rev_p)
		diff=abs(fwd_p-rev_p)
		if pvalue<=cutoff:
			r.append((l[0],l[1],str(pvalue),str(diff)))
		#print([fwd_m,rev_m,fwd_um,rev_um],fwd_p,rev_p,pvalue)
	return r
if __name__ == '__main__':
	prefix=sys.argv[1]
	pvalues=float(sys.argv[2])
	f=sys.argv[3] #CX_report.txt.gz list
	global epsilon
	epsilon=1e-07
	Fr_CG=open(prefix+"_CG_hemi_ori.bed","w")
	Fr_CHG=open(prefix+"_CHG_hemi_ori.bed","w")
	rCG=GetCGHemi(f,pvalues)
	rCHG=GetCHGHemi(f,pvalues)
	for ri in rCG:
		Fr_CG.write(ri[0]+"\t"+ri[1]+"\t"+str(int(ri[1])+1)+"\t.\t"+str(ri[2])+"\t"+str(ri[3])+"\n")
	for ri in rCHG:
		Fr_CHG.write(ri[0]+"\t"+ri[1]+"\t"+str(int(ri[1])+1)+"\t.\t"+str(ri[2])+"\t"+str(ri[3])+"\n")
	Fr_CG.close()
	Fr_CHG.close()
	os.system("bedtools sort -i %s_CG_hemi_ori.bed|bedtools merge -d 5 -i - >%s_CG_hemi_final.bed"%(prefix,prefix)) #-d need change
	os.system("bedtools sort -i %s_CHG_hemi_ori.bed|bedtools merge -d 5 -i - >%s_CHG_hemi_final.bed"%(prefix,prefix)) #-d need change
