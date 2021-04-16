#!/usr/bin/env python
import sys,pandas,os,collections
import numpy as np
import statsmodels.api as sm
import scipy.stats as stats
def FisherExactTest(MList):
	pass
def Gmt(MList):
	X=[]
	y=[]
	epsilon=1e-07
	if len(set(MList))==1:
		return 1
	#print(MList)
	MList=np.array(MList)
	MList=np.clip(MList,epsilon,1-epsilon)
	y=np.log(MList/(1-MList))
	#print(y)
	X=[1 if i%2==0 else 0 for i in range(1,len(y)+1)]
	#print(X)
	X=sm.add_constant(X,prepend=False)
	mod=sm.OLS(y,X)
	res = mod.fit()
	#print(res.summary())
	p=res.pvalues[0]
	#print(p)
	return p
def GetCGHemi(f,d,dr):
	lines=os.popen("zgrep -P '\tCG\t' %s"%f).readlines()
	for i,x in enumerate(lines):
		x=x.rstrip()
		l=x.split("\t")
		if not l[0] in d:
			d[l[0]]={}
		if not l[0] in dr:
			dr[l[0]]={}
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
			dr[l[0]][l[1]].append(fwd_m+fwd_um)
			if rev_m+rev_um!=0:
				d[l[0]][l[1]].append(rev_m/(rev_m+rev_um))
			else:
				d[l[0]][l[1]].append(0)
			dr[l[0]][l[1]].append(rev_m+rev_um)
		except:
			d[l[0]][l[1]]=[]
			dr[l[0]][l[1]]=[]
			if fwd_m+fwd_um!=0:
				d[l[0]][l[1]].append(fwd_m/(fwd_m+fwd_um))
			else:
				d[l[0]][l[1]].append(0)
			dr[l[0]][l[1]].append(fwd_m+fwd_um)
			if rev_m+rev_um!=0:
				d[l[0]][l[1]].append(rev_m/(rev_m+rev_um))
			else:
				d[l[0]][l[1]].append(0)
			dr[l[0]][l[1]].append(rev_m+rev_um)
	return d,dr
def GetCHGHemi(f,d,dr):
	lines=os.popen("zgrep -P '\tCHG\t' %s"%f).readlines()
	for i,x in enumerate(lines):
		x=x.rstrip()
		l=x.split("\t")
		if not l[0] in d:
			d[l[0]]={}
		if not l[0] in dr:
			dr[l[0]]={}
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
			dr[l[0]][l[1]].append(fwd_m+fwd_um)
			if rev_m+rev_um!=0:
				d[l[0]][l[1]].append(rev_m/(rev_m+rev_um))
			else:
				d[l[0]][l[1]].append(0)
			dr[l[0]][l[1]].append(rev_m+rev_um)
		except:
			d[l[0]][l[1]]=[]
			dr[l[0]][l[1]]=[]
			if fwd_m+fwd_um!=0:
				d[l[0]][l[1]].append(fwd_m/(fwd_m+fwd_um))
			else:
				d[l[0]][l[1]].append(0)
			dr[l[0]][l[1]].append(fwd_m+fwd_um)
			if rev_m+rev_um!=0:
				d[l[0]][l[1]].append(rev_m/(rev_m+rev_um))
			else:
				d[l[0]][l[1]].append(0)
			dr[l[0]][l[1]].append(rev_m+rev_um)
	return d,dr
if __name__ == '__main__':
	prefix=sys.argv[1]
	pvalues=float(sys.argv[2])
	flist=sys.argv[3:] #CX_report.txt.gz list
	diff=0.75 #need change
	cutoff=4 #need change
	Fr_CG=open(prefix+"_CG_hemi_ori.bed","w")
	Fr_CHG=open(prefix+"_CHG_hemi_ori.bed","w")
	dCHG={}
	drCHG={}
	dCG={}
	drCG={}
	for f in flist:
		dCG,drCG=GetCGHemi(f,dCG,drCG)
		dCHG,drCHG=GetCHGHemi(f,dCHG,drCHG)
	for c in dCG:
		for p in dCG[c]:
			MList=dCG[c][p]
			#diff=np.mean([abs(MList[0]-MList[1]),abs(MList[2]-MList[3])])
			fwd=np.mean([MList[0],MList[2]])
			rev=np.mean([MList[1],MList[3]])
			p_value=Gmt(MList)
			#print(c,p,dCG[c][p],p_value)
			ReadsL=drCG[c][p]
			fwdRs=(ReadsL[0]+ReadsL[2])/2.
			revRs=(ReadsL[1]+ReadsL[3])/2.
			if p_value<=pvalues:
				Fr_CG.write(c+"\t"+p+"\t"+str(int(p)+1)+"\t"+str(fwd)+"\t"+str(rev)+"\t"+str(fwdRs)+"\t"+str(revRs)+"\n")
	for c in dCHG:
		for p in dCHG[c]:
			MList=dCHG[c][p]
			fwd=np.mean([MList[0],MList[2]])
			rev=np.mean([MList[1],MList[3]])
			p_value=Gmt(MList)
			ReadsL=drCHG[c][p]
			fwdRs=(ReadsL[0]+ReadsL[2])/2.
			revRs=(ReadsL[1]+ReadsL[3])/2.
			if p_value<=pvalues:
				Fr_CHG.write(c+"\t"+p+"\t"+str(int(p)+2)+"\t"+str(fwd)+"\t"+str(rev)+"\t"+str(fwdRs)+"\t"+str(revRs)+"\n")
	Fr_CG.close()
	Fr_CHG.close()
	os.system("awk '($4-$5>=%s||$4-$5<=%s)&&($6>=%s||$7>=%s){print $1\"\t\"$2\"\t\"$3}' %s_CG_hemi_ori.bed|bedtools sort -i -|bedtools merge -d 5 -i - >%s_CG_hemi_final.bed"%(diff,-1*diff,cutoff,cutoff,prefix,prefix)) #-d need change
	os.system("awk '($4-$5>=%s||$4-$5<=%s)&&($6>=%s||$7>=%s){print $1\"\t\"$2\"\t\"$3}' %s_CHG_hemi_ori.bed|bedtools sort -i -|bedtools merge -d 5 -i - >%s_CHG_hemi_final.bed"%(diff,-1*diff,cutoff,cutoff,prefix,prefix)) #-d need change
