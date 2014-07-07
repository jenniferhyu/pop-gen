import numpy
import scipy
import scipy.stats as stat
import random
from math import sqrt
import sys
import decimal
#############################
#Make R equivalent commands #
#############################
rnorm = numpy.random.normal
pnorm = stat.norm.cdf
sd = numpy.std
mean = numpy.mean
chisq = stat.chisquare
############################
fileStart=str(sys.argv[1])
replicates=int(sys.argv[2])
cutoff=1-((float(sys.argv[3]))/100)
control=bool(sys.argv[4]) #only when lambda=0 is this True.
roundoff=decimal.Decimal(str(cutoff))
roundoff=-(roundoff.as_tuple().exponent)
rep=0
logged=[]
logged250=[]
logged500=[]
pheno=numpy.genfromtxt("phenotypes."+fileStart,delimiter="\t")
Vg=pheno[:,0]
Ve=pheno[:,1]
Vt=numpy.add(Vg,Ve)
if not control:
	PV=numpy.genfromtxt("effectfiletest."+fileStart,dtype=float,delimiter="\t")
	#positions=PV[:,0]
	values=numpy.matrix(PV[:,2])
	values=values.T
	causes=numpy.genfromtxt("matrixprep"+fileStart+".txt",delimiter=" ")
	haplo=numpy.dot(causes,values)
else:
	haplo=numpy.matrix([0]*40000)
	haplo=haplo.T
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
		ppick=list(set(pick)-set(mpick)) #randomized
		maternal=[]
		paternal=[]
		for m in mpick:
			maternal.append(parents[m])
		for p in ppick:
			paternal.append(parents[p])	
		children=numpy.empty(shape=(2,2))
		for i in range(2):
			for j in range(2):
				children[i,j]=sqrt(maternal[i]*paternal[j])
		children=numpy.asarray(children).reshape(-1)
		select=[random.randint(0,3) for select in range(2)]
		siblings=[]
		for i in select:
			siblings.append(children[i])
		noise=rnorm(scale=sd(Ve),size=2)
		siblings=numpy.add(siblings,noise)
		zscores=[]
		k=0
		while(k<len(siblings)):		
			score=round(pnorm(siblings[k],loc=mean(Vt),scale=sd(Vt)),roundoff)
			if(score>=cutoff): 
				zscores.append(select[k])
			k+=1
		if len(zscores)==2:
			ibd=0
			if(zscores[0]==zscores[1]):
				N2+=1
			elif (sum(zscores)==3):
				N0+=1
			else:
				N1+=1
			zscores=[] #sanity reset
			Ntotal+=1
			observed=numpy.array([N0,N1,N2])
			expected=numpy.array([250,500,250])
			if Ntotal==250:
				pvalue250=chisq(observed,expected)[1]
				logged250.append(pvalue250)
			if Ntotal==500:
				pvalue500=chisq(observed,expected)[1]
				logged500.append(pvalue500)
	IBD=numpy.array([N0,N1,N2])
	print(IBD)
	pvalue=chisq(IBD,expected)[1]
	logged.append(pvalue)
	rep+=1
N=[[250]*len(logged250)]
N.append([500]*len(logged500))
N.append([1000]*len(logged))
N=[item for sublist in N for item in sublist]
allLogs=[logged250,logged500,logged]
allLogs=[item for sublist in allLogs for item in sublist]
allLogs=-numpy.log10(allLogs)
print(log)
