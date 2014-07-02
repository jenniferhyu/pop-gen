import numpy
import scipy
import scipy.stats as stat
#import argparse
from random import randint
from math import sqrt
# parser = argparse.ArgumentParser(description="Paramters")
# parser.add_argument('file#',metavar="f",type=str,dest="opened_file")
# parser.add_argument('replicates',metavar='replicates',type=int)
# parser.add_argument('cutoff', type=float)
#############################
# Make R-equivalent commands#
#############################
rnorm = stat.norm.rvs
pnorm = stat.norm.cdf
sd = stat.std
mean = stat.mean
chisq = stat.chisquare
############################
fileStart=str(sys.argv[1])
replicates=int(sys.argv[2])
cutoff=1-((float(sys.argv[3]))/100)
r=1
logged=[]
logged250=[]
logged500[]
PV=numpy.genfromtxt("effectfiletest."+fileStart,delimiter="\t")
positions=PV[:,0]
values=PV[:,1]
pheno=numpy.genfromtxt("phenotypes."+fileStart,delimter="\t")
Vg=PV[:,0]
Ve=PV[:,1]
Vt=numpy.sum(Vg,Ve,axis=0)
causes=numpy.loadtxt("matrixprep"+fileStart+".txt",delimiter=" ")
haplo=causes.dot(values)
while(r<=replicates):
	Ntotal=0
	NTotalGoal=1000
	N0=0
	N1=0
	N2=0
	while(Ntotal<NTotalGoal):
		gametes=[randint(1,4000) for p in range(0,4)]
		parents=[]
		for g in gametes:
			parents.append(haplo[g])
		pick=[randint(1,4) for p in range(0,4)]
		mpick=pick[0:2]
		ppick=pick[2:4]
		maternal=[]
		paternal=[]
		for m in mpick:
			maternal.append(parents[m])
		for p in ppick:
			paternal.append(parents[p])
		children=numpy.empty(shape=(2,2))
		i=0
		j=0
		while(i<2):
			while(j<2):
				children[i,j]=sqrt(maternal[i]*paternal[j])
				j+=1
			i+=1
		select=[randint(1,2) for select in range(0,4)]
		siblings=[]
		k=0
		while(k+1<len(select)):
			siblings.append(children[siblings[k],siblings[k+1]])
			k+=2
		noise=list(rnorm(loc=0,scale=sd(Ve),size=4))
		siblings=numpy.add(siblings,noise)
		zscores=[]
		for s in siblings:
			if(round(pnorm(s,mean(Vt),sd(Vt)),3)>=cutoff):
				zscores.append(s)
		if len(zscores)==2:
			alleles=numpy.reshape(select,(2,-1)) #-1????
			ibd=0
			for r in alleles:
				if r[0]==r[1]:
					ibd+=2
			if (ibd==0):
				N0+=1
			if (ibd==1):
				N1+=1
			if (ibd==2):
				N2+=1
		zscores=[]
		Ntotal+=1
		observed=numpy.array([N0.,N1.,N2.])
		expected=numpy.array([.25,.50,.25])*numpy.sum(observed)
		if Ntotal==250:
			pvalue250=chisq(observed,expected)[1]
			logged250.append(pvalue250)
		if Ntotal=500:
			pvalue500=chisq(observed,expected)[1]
			logged500.append(pvalue500)
	IBD=numpy.array([N0.,N1.,N2.])
	pvalue=chisq(observed,expected)[1]
	logged.append(pvalue)
	r+=1
	

def ibd_count(x):
	if x[,1] == x[,2]:
		return 1
	else:
		return 0




