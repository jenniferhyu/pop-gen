import numpy
import scipy
import scipy.stats as stat
import random
from math import sqrt
import sys
import decimal
import itertools
#############################
#Make R equivalent commands #
#############################
rnorm = stat.norm.rvs
pnorm = stat.norm.cdf
sd = numpy.std
mean = numpy.mean
chisq = stat.chisquare
############################
fileStart=str(sys.argv[1])
replicates=int(sys.argv[2])
cutoff=1-((float(sys.argv[3]))/100)
roundoff=decimal.Decimal(str(cutoff))
roundoff=-(roundoff.as_tuple().exponent)
rep=0
logged=[]
logged250=[]
logged500=[]
PV=numpy.genfromtxt("effectfiletest."+fileStart,dtype=float,delimiter="\t")
#positions=PV[:,0]
values=numpy.matrix(PV[:,2])
values=values.transpose()
pheno=numpy.genfromtxt("phenotypes."+fileStart,delimiter="\t")
Vg=pheno[:,0]
Ve=pheno[:,1]
Vt=numpy.add(Vg,Ve)
causes=numpy.genfromtxt("matrixprep"+fileStart+".txt",delimiter=" ")
haplo=numpy.dot(causes,values)
while(rep<replicates):
	Ntotal=0
	NTotalGoal=1000
	N0=0
	N1=0
	N2=0
	while(Ntotal<NTotalGoal):
		gametes=random.sample(range(40000),4)
		parents=[]
		for g in gametes:
			parents.append(haplo[g,0])
		pick=random.sample(range(4),4)
		mpick=random.sample(range(4),2)
		ppick=list(set(pick)-set(mpick))
		maternal=[]
		paternal=[]
		for m in mpick:
			maternal.append(parents[m])
		for p in ppick:
			paternal.append(parents[p])	
		children=numpy.empty(shape=(2,2))
		for i in range(0,2):
			for j in range(0,2):
				children[i,j]=sqrt(maternal[i]*paternal[j])
		select=[random.randint(0,1) for select in range(0,4)]
		siblings=[]
		k=0
		while(k<len(select)):
			siblings.append(children.item(select[k],select[k+1]))
			k+=2
		noise=numpy.random.normal(scale=sd(Ve),size=2)
		siblings=numpy.add(siblings,noise)
		zscores=[]
		for s in siblings:
			score=round(pnorm(s,loc=mean(Vt),scale=sd(Vt)),roundoff)
			if(score>=cutoff): 
				zscores.append(s)
		if len(zscores)==2:
			alleles=numpy.reshape(select,(2,-1))
			ibd=0
			for i in range(0,2):
				if(alleles[0,i]==alleles[1,i]):
					ibd+=1
			if (ibd==0):
				N0+=1
			if (ibd==1):
				N1+=1
			if (ibd==2):
				N2+=1
			zscores=[]
			Ntotal+=1
			observed=numpy.array([N0,N1,N2])
			expected=numpy.array([250,500,250])
			#expected=numpy.array([250,500,250])
			if Ntotal==250:
				pvalue250=chisq(observed,expected)[1]
				logged250.append(pvalue250)
			if Ntotal==500:
				pvalue500=chisq(observed,expected)[1]
				logged500.append(pvalue500)
	IBD=numpy.array([N0,N1,N2])
	pvalue=chisq(IBD,expected)[1]
	logged.append(pvalue)
	rep+=1
N=[[250]*len(logged250)]
N.append([500]*len(logged500))
N.append([1000]*len(logged))
N_flat=[item for sublist in N for item in sublist]
log=-numpy.log10(logged)
print(log)