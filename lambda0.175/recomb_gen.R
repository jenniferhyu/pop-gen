#!/usr/bin/Rscript
decimalPlaces<-function(x){
	if ((x %% 1) != 0) {
        nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
    } else {
        return(0)
    }
}
# start.time<-Sys.time()
#args<-commandArgs(trailingOnly=TRUE)
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
rep<-1
replicate<-1
places<-21 #should be an odd number
totalMarkers<-20
haploLength<-100
sp<-haploLength/(places-1) 
pos<-seq(0,haploLength,sp)
#pos<-pos[-((length(pos)%/%2)+1)]
#space<-(haploLength/(places-1))/100
roundcM<-function(n, center=c, upper=up, spacing=space){
	if (n<upper && n>=center) {
		return((places%/%2)+1) #round down
	} else if (n>=upper) {
		return(floor(n/spacing)+2)
	} else {
		return(floor(n/spacing)+1)
	}
}
other<-function(tb){
	return(x[!x%in%tb])
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
result<-vector()
famTable<-vector()	
while(rep<=replicate){
	Ntotal<-0
	#NtotalGoal<-as.numeric(args[1])
	NtotalGoal<-250 #possibly moved to command-line argument later
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
		parentScores<-vector()
		i=1
		while(i+1<=length(select)){			
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
			half<-places%/%2
			msample1<-c(sample(1:5,totalMarkers,T))
			msample2<-c(sample(1:5,totalMarkers,T))
			M1<-double(places)
			M1[1:half]<-msample1[1:half]
			M1[(ceiling(places/2)+1):places]<-msample1[(half+1):length(msample1)]
			M1[ceiling(places/2)]<-maternal[1]
			M1<-M1[1:21]
			M2<-double(places)
			M2[1:half]<-msample2[1:half]
			M2[(ceiling(places/2)+1):places]<-msample2[(half+1):length(msample2)]
			M2[ceiling(places/2)]<-maternal[2]
			M2<-M2[1:21]
			M<-rbind(M1,M2)
			M<-M[,colSums(is.na(M))==0]
			psample1<-c(sample(1:5,totalMarkers,T))
			psample2<-c(sample(1:5,totalMarkers,T))
			P1<-double(places)
			P1[1:half]<-psample1[1:half]
			P1[(ceiling(places/2)+1):places]<-psample1[(half+1):length(psample1)]
			P1[ceiling(places/2)]<-paternal[1]
			P1<-P1[!is.na(P1)]
			P2<-double(places)
			P2[1:half]<-psample2[1:half]
			P2[(ceiling(places/2)+1):places]<-psample2[(half+1):length(psample2)]
			P2[ceiling(places/2)]<-paternal[2]
			P2<-P2[!is.na(P2)]		
			P<-rbind(P1,P2)
			P<-P[,colSums(is.na(P))==0]
			space<-(haploLength/(places-1))/100 #or haploLength/places
			c<-haploLength/2
			up<-pos[match(c,pos)+1]
			site<-(length(M1)%/%2)+1
			#################
			# Recombination #
			#################
			counter<-1
			founder<-rbind(M,P)
			pipe<-vector()
			set<-vector()
			while(counter<=2) {
				recombM<-FALSE
				recombP<-FALSE
				rsite1<-runif(1)
				rsite2<-runif(1)
				tb1<-sample(1:2,1) #start from top or bottom
				tb2<-sample(1:2,1)
				if(counter==1){
					pair<-select[1:2]
				} else {
					pair<-select[3:4]
				}
				if (rsite1<=haploLength/100) {
					recombM<-TRUE
					cut1<-roundcM(rsite1)
					mom1<-c(M[tb1,][1:cut1],M[other(tb1),][(cut1+1):length(M[1,])])
					mom2<-c(M[other(tb1),][1:cut1],M[tb1,][(cut1+1):length(M[1,])])
					rM<-rbind(mom1,mom2)				
					rM<-rM[,colSums(is.na(rM))==0]
					splitM<-rbind(c(tb1,other(tb1)),c(other(tb1),tb1))
				}				
				if (rsite2<=haploLength/100){
					recombP<-TRUE
					cut2<-roundcM(rsite2)				
					dad1<-c(P[tb2,][1:cut2],P[other(tb2),][(cut2+1):length(P[1,])])
					dad2<-c(P[other(tb2),][1:cut2],P[tb2,][(cut2+1):length(P[1,])])
					rP<-rbind(dad1,dad2)
					rP<-rP[,colSums(is.na(rP))==0]
					splitP<-rbind(c(tb2,other(tb2)),c(other(tb2),tb2))
				}
				if (recombP&&recombM){
					person<-rbind(rM[pair[1],],rP[pair[2],])
					pinfo<-c(splitM[pair[1],],splitP[pair[2],])
				} else if (recombM){
					person<-rbind(rM[pair[1],],P[pair[2],])
					pinfo<-c(splitM[pair[1],],splitP[pair[2],])
				} else if (recombP) {
					person<-rbind(M[pair[1],],rP[pair[2],])
					pinfo<-c(splitM[pair[1],],splitP[pair[2],])
				} else {
					person<-rbind(M[pair[1],],P[pair[2],])	
					pinfo<-c(splitM[pair[1],],splitP[pair[2],])				
				}
				#person<-rbind(rM[pair[1],],rP[pair[2],])
				rownames(person)<-c(paste("Sib",counter,"-Mom",sep=""),paste("Sib",counter,"-Dad",sep=""))
				string<-c(pinfo[1],rsite1,pinfo[2],pinfo[3],rsite2,pinfo[4],noise[counter])
				pipe<-append(pipe,string)
				set<-rbind(set,person)
				counter<-counter+1
			}
			result<-rbind(founder,set)
			output<-rbind(output,result)
			#output<-cbind(output,occurrence)
			zscores<-vector()
			info<-c(maternal,paternal,
				pipe)
			famTable<-rbind(famTable,info)
		}			
	}
	rep<-rep+1
}
IBD<-c(N0,N1,N2)
p.value<-chisq.test(IBD,p=c(1,2,1),rescale.p=TRUE)$p.value
LODscore<-(-log10(p.value))
sink("hidden_LOD20.txt",append=TRUE)
print(LODscore)
sink()
write.table(famTable,file="family_info.fam",row.names=FALSE,col.names=FALSE)
# recomb_check<-sum(occurrence==TRUE)/length(occurrence)
# print(recomb_check)
write.table(output,file="foundertest2.txt",row.names=TRUE,col.names=FALSE)
write.table(output,file="no_row2.txt",row.names=FALSE,col.names=FALSE)
