#!/usr/bin/Rscript
decimalPlaces<-function(x){
	if ((x %% 1) != 0) {
        nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
    } else {
        return(0)
    }
}
# start.time<-Sys.time()
# args<-commandArgs(trailingOnly=TRUE)
# if (length(args)<3) {
# 	warning("Must have 3 args: file, replicate, and cutoff")
# }
# fileStart<-as.numeric(args[1])
# replicate<-as.numeric(args[2])
# cutoff<-(1-(as.numeric(args[3])/100))
# places<-as.numeric(args[4])
# haploLength<-as.numeric(args[5])
# r<-1
# roundoff<-decimalPlaces(cutoff)
# ###################
# Manual switches  #
##################
cutoff<-1-(7.5/100)
roundoff<-decimalPlaces(cutoff)
fileStart<-1
r<-1
replicate<-1
places<-21
haploLength<-100
sp<-haploLength/(places-1)
pos<-seq(0,haploLength,sp)
roundcM<-function(n, center=c, upper=up, spacing=space){
	if (n<upper && n>=center) {
		return((places%/%2)+1)
	} else if (n>=upper) {
		return(floor(n/spacing)+2)
	} else {
		return(floor(n/spacing)+1)
	}
}
###########
# Loading #
###########
PV <- read.table(paste("effectfiletest.",as.character(fileStart),sep=""), header= FALSE, sep="\t")
positions<-PV$V1
values<-PV$V3
pheno<-read.table(paste("phenotypes.",as.character(fileStart),sep=""), header=FALSE, sep="\t")
Vg<-pheno$V1
Ve<-pheno$V2
Vt<-Vg+Ve	
causes<-as.matrix(read.table(paste("matrixprep",as.character(fileStart),".txt",sep=""),header=FALSE, sep=" "))
causes<-causes[,colSums(is.na(causes))==0] #removes any NA values that may arise
colnames(causes)<-positions
haplo<-causes%*%values
x<-c(1,2)
y<-c(1,2)
output<-vector()	
while(r<=replicate){
	Ntotal<-0
	NtotalGoal<-1000 #possibly moved to command-line argument later
	N0<-0
	N1<-0
	N2<-0
	#print(paste("Read file",as.character(fileStart),"on replicate",as.character(r)))		
	while(Ntotal<NtotalGoal){				
		gametes<-sample(1:40000, 4)
		parents <- vector()
		for (g in gametes) {			
			parents<-append(parents,haplo[g])
		}
		pick<-sample(1:4,4)
		mpick<-sample(1:4,2)
		ppick<-pick[!pick%in%mpick]
		maternal<-vector()
		paternal<-vector()
		#cM<-(c(0,9,18,22.5,27,36,45)/100)
		for (m in mpick) {			
			maternal<-append(maternal,parents[m])
		}
		for (p in ppick){			
			paternal<-append(paternal,parents[p])
		}
		select<-sample(1:2,4,replace=TRUE)
		siblings<-vector() #randomly choosing two siblings with only Vg
		i=1
		while(i+1<=length(select)){			
			#siblings<-append(siblings,sqrt(M[select[i],][4]*P[select[i+1],][4]))
			siblings<-append(siblings,sqrt(maternal[select[i]]*paternal[select[i+1]]))
			i<-i+2
		}
		noise <- rnorm(2, 0, sd(Ve)) #random number generation
		siblings<-siblings+noise		
		zscores<-ifelse(round(pnorm(siblings,mean(Vt),sd(Vt)),roundoff)>=cutoff, TRUE, FALSE)
		if (sum(zscores==TRUE)==2){
			Ntotal<-Ntotal+1
			if (select[1]==select[3]&&select[2]==select[4]){
				N2<-N2+1
			} else if (select[1]==select[3]||select[2]==select[4]){
				N1<-N1+1
			} else {
				N0<-N0+1
			}
			msample1<-c(sample(1:5,places-1,T))
			msample2<-c(sample(1:5,places-1,T))
			M1<-double(places)
			M1[1:places%/%2]<-msample1[1:places%/%2]
			M1[ceiling(places/2)+1:places]<-msample1[(places%/%2)+1:places-1]
			M1[ceiling(places/2)]<-maternal[1]
			M2<-double(places)
			M2[1:places%/%2]<-msample2[1:places%/%2]
			M2[ceiling(places/2)+1:places]<-msample2[(places%/%2)+1:places-1]
			M2[ceiling(places/2)]<-maternal[2]			
			#M1<-c(sample(1:5,places%/%2,replace=TRUE),maternal[1],sample(1:5,places%/%2,replace=TRUE))
			#M2<-c(sample(1:5,places%/%2,replace=TRUE),maternal[2],sample(1:5,places%/%2,replace=TRUE))
			M<-rbind(M1,M2)
			M<-M[,colSums(is.na(M))==0]
			psample1<-c(sample(1:5,places-1,T))
			psample2<-c(sample(1:5,places-1,T))
			P1<-double(places)
			P1[1:places%/%2]<-psample1[1:places%/%2]
			P1[ceiling(places/2)+1:places]<-psample1[(places%/%2)+1:places-1]
			P1[ceiling(places/2)]<-paternal[1]
			P2<-double(places)
			P2[1:places%/%2]<-psample2[1:places%/%2]
			P2[ceiling(places/2)+1:places]<-psample2[(places%/%2)+1:places-1]
			P2[ceiling(places/2)]<-paternal[2]		
			#P1<-c(sample(1:5,places%/%2,replace=TRUE),paternal[1],sample(1:5,places%/%2,replace=TRUE))
			#P2<-c(sample(1:5,places%/%2,replace=TRUE),paternal[2],sample(1:5,places%/%2,replace=TRUE))
			P<-rbind(P1,P2)
			P<-P[,colSums(is.na(P))==0]
			space<-haploLength/(places-1) #or haploLength/places
			c<-haploLength/2
			up<-pos[match(c,pos)+1]
			site<-(length(M1)%/%2)+1
			#################
			# Recombination #
			#################
			recombM<-FALSE
			recombP<-FALSE
			rsite1<-runif(1)
			rsite2<-runif(1)
			tb1<-sample(1:2,1) #start from top or bottom
			tb2<-sample(1:2,1)
			other1<-x[!x%in%tb1]
			other2<-y[!y%in%tb2]
			if (rsite1<=haploLength/100) {
				recombM<-TRUE
				cut1<-roundcM(rsite1)
				mom1<-c(M[tb1,][1:cut1],M[other1,][cut1+1:length(M1)])
				mom2<-c(M[other1,][1:cut1],M[tb1,][cut1+1:length(M1)])
				rM<-rbind(mom1,mom2)				
				rM<-rM[,colSums(is.na(rM))==0]
			}
			if (rsite2<=haploLength/100){
				recombP<-TRUE
				cut2<-roundcM(rsite2)				
				dad1<-c(P[tb2,][1:cut2],P[other2,][cut2+1:length(P1)])
				dad2<-c(P[other2,][1:cut2],P[tb2,][cut2+1:length(P1)])
				rP<-rbind(dad1,dad2)
				rP<-rP[,colSums(is.na(rP))==0]
			}
			founder<-rbind(M,P)
			#rsites<-append(rsites, rsite1)
			#rsites<-append(rsites, rsite2)
			if (recombM){
				person1<-rbind(rM[select[1],],P[select[2],])
			} else if (recombP) {
				person1<-rbind(M[select[1],],rP[select[2],])
			} else if (recombP&&recombM) {
				person1<-rbind(rM[select[1],],rP[select[2],])
			} else {
				person1<-rbind(M[select[1],],P[select[2],])
			}
			#person1[,site]<-siblings[1]
			rownames(person1)<-c("Sib1-Mom","Sib1-Dad")
			if (recombM){
				person2<-rbind(rM[select[3],],P[select[4],])
			} else if (recombP) {
				person2<-rbind(M[select[3],],rP[select[4],])
			} else if (recombP&&recombM) {
				person2<-rbind(rM[select[3],],rP[select[4],])
			} else {
				person2<-rbind(M[select[3],],P[select[4],])
			}
			#person2[,site]<-siblings[2]
			rownames(person2)<-c("Sib2-Mom","Sib2-Dad")
			result<-rbind(founder,person1,person2)
			output<-rbind(output,result)
			#output<-cbind(output,occurrence)
			zscores<-vector()
		}			
	}
	r<-r+1
}
IBD<-c(N0,N1,N2)
p.value<-chisq.test(IBD,p=c(1,2,1),rescale.p=TRUE)$p.value
LODscore<-(-log10(p.value))
print(LODscore)
if(LODscore>=3){
sink("hidden_LOD20.txt",append=TRUE)
print(LODscore)
sink()
}
# recomb_check<-sum(occurrence==TRUE)/length(occurrence)
# print(recomb_check)
#sink("recomb_rates.txt",append=TRUE)
#print(recomb_check)
#sink()
write.table(output,file="foundertest2.txt",row.names=TRUE,col.names=FALSE)
write.table(output,file="no_row2.txt",row.names=FALSE,col.names=FALSE)
#write.table(rsites, file="rsites.txt",col.names=FALSE)