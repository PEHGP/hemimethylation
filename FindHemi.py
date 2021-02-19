#!/usr/bin/env python
import sys,pandas,os,collections
import numpy as np
import statsmodels.api as sm
def FisherExactTest():
	pass
def Gmt(MList,treatment=1):
	X=[]
	y=[]
	if len(set(MList))==1:
		return 1
	for m in MList:
		if m==1:
			m=0.99999999
		elif m==0:
			m=0.00000001
		y.append(np.log(m/(1-m)))
		X.append(treatment)
	X=sm.add_constant(X,prepend=False)
	mod=sm.OLS(y,X)
	res = mod.fit()
	#print(res.summary())
	p=res.pvalues[0]
	return p
def GetCGHemi(f,d):
	lines=os.popen("zgrep -P '\tCG\t' %s"%f).readlines()
	for i,x in enumerate(lines):
		x=x.rstrip()
		l=x.split("\t")
		if not l[0] in d:
			d[l[0]]={}
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
		try:
			if fwd_m+fwd_um!=0:
				d[l[0]][l[1]].append(fwd_m/(fwd_m+fwd_um))
			else:
				d[l[0]][l[1]].append(0)
			if rev_m+rev_um!=0:
				d[l[0]][l[1]].append(rev_m/(rev_m+rev_um))
			else:
				d[l[0]][l[1]].append(0)
		except:
			d[l[0]][l[1]]=[]
			if fwd_m+fwd_um!=0:
				d[l[0]][l[1]].append(fwd_m/(fwd_m+fwd_um))
			else:
				d[l[0]][l[1]].append(0)
			if rev_m+rev_um!=0:
				d[l[0]][l[1]].append(rev_m/(rev_m+rev_um))
			else:
				d[l[0]][l[1]].append(0)
	return d
def GetCHGHemi(f,d):
	lines=os.popen("zgrep -P '\tCHG\t' %s"%f).readlines()
	for i,x in enumerate(lines):
		x=x.rstrip()
		l=x.split("\t")
		if not l[0] in d:
			d[l[0]]={}
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
		try:
			if fwd_m+fwd_um!=0:
				d[l[0]][l[1]].append(fwd_m/(fwd_m+fwd_um))
			else:
				d[l[0]][l[1]].append(0)
			if rev_m+rev_um!=0:
				d[l[0]][l[1]].append(rev_m/(rev_m+rev_um))
			else:
				d[l[0]][l[1]].append(0)
		except:
			d[l[0]][l[1]]=[]
			if fwd_m+fwd_um!=0:
				d[l[0]][l[1]].append(fwd_m/(fwd_m+fwd_um))
			else:
				d[l[0]][l[1]].append(0)
			if rev_m+rev_um!=0:
				d[l[0]][l[1]].append(rev_m/(rev_m+rev_um))
			else:
				d[l[0]][l[1]].append(0)
	return d
if __name__ == '__main__':
	prefix=sys.argv[1]
	pvalues=float(sys.argv[2])
	flist=sys.argv[3:] #CX_report.txt.gz list
	Fr_CG=open(prefix+"_CG_hemi_ori.bed","w")
	Fr_CHG=open(prefix+"_CHG_hemi_ori.bed","w")
	dCHG={}
	dCG={}
	for f in flist:
		dCG=GetCGHemi(f,dCG)
		dCHG=GetCHGHemi(f,dCHG)
	for c in dCG:
		for p in dCG[c]:
			 p_value=Gmt(dCG[c][p])
			 if p_value<=pvalues:
			 	Fr_CG.write(c+"\t"+p+"\t"+str(int(p)+1)+"\n")
	for c in dCHG:
		for p in dCHG[c]:
			p_value=Gmt(dCHG[c][p])
			if p_value<=pvalues:
				Fr_CHG.write(c+"\t"+p+"\t"+str(int(p)+2)+"\n")
	Fr_CG.close()
	Fr_CHG.close()
	os.system("bedtools sort -i %s_CG_hemi_ori.bed|bedtools merge -d 5 -i - >%s_CG_hemi_final.bed"%(prefix,prefix)) #-d need change
	os.system("bedtools sort -i %s_CHG_hemi_ori.bed|bedtools merge -d 5 -i - >%s_CHG_hemi_final.bed"%(prefix,prefix)) #-d need change
