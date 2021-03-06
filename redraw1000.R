#!/usr/bin/Rscript
args<-commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
	warning("Must have 3 args: file, replicate, and cutoff")
}
fileStart<-as.numeric(args[1])
replicate<-as.numeric(args[2])
cutoff<-(1-(as.numeric(args[3])/100))
template<-paste("result-",args[3],"-",fileStart,".txt",sep="")
r<-1
logged<-vector()
logged250<-vector()
logged500<-vector()
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
	while(r<=replicate){
		Ntotal<-0
		NtotalGoal<-1000
		N0<-0
		N1<-0
		N2<-0
		#^every new replicate, I need to clear the pre-existing N-associated values^#
		#print(paste("Read file",as.character(fileStart),"on replicate",as.character(r)))		
		while(Ntotal<NtotalGoal){				
			#############################
			#to generate IBDs ultimately#
			#############################			
			gametes<-sample(1:40000, 4)
			parents <- vector()
			for (i in 1:4) {			
				parents<-append(parents,haplo[gametes[i]])
			}
			pick<-sample(1:4,4)
			mpick<-sample(1:4,2)
			ppick<-pick[!pick%in%mpick]
			maternal<-vector()
			paternal<-vector()				
			for (i in 1:2) {			
				maternal<-append(maternal,parents[mpick[i]])
			}
			for (i in 1:2){			
				paternal<-append(paternal,parents[ppick[i]])
			}
			children<-matrix(NA,2,2)
			for (i in 1:2) {
				for (j in 1:2) {
					children[i,j]<-sqrt(maternal[i]*paternal[j])
				}
			}			
			select<-sample(1:2,4,replace=TRUE)	
			siblings<-vector() #randomly choosing two siblings with only Vg
			i=1
			while(i+1<=length(select)){			
				siblings<-append(siblings,children[select[i],select[i+1]])
				i<-i+2
			}
			noise <- rnorm(2, 0, sd(Ve)) #random number generation
			siblings<-siblings+noise		
			zscores<-vector()
			for (i in 1:length(siblings)){
				if(round(pnorm(siblings[i],mean(Vt),sd(Vt)),2)>=cutoff){ #top 15% of phenotypic variance (previously hard-coded)				
					zscores<-append(zscores,siblings[i])
				}
			}
			if (length(zscores)!=2){
				rm(children,parents,gametes,maternal,paternal,siblings,pick)										
			} else {									
				alleles<-matrix(select,2,2)				
				ibd = 0 #resetting ibd for each iteration of the loop
				for (i in 1:nrow(alleles)) {
					if (alleles[i,1]==alleles[i,2]){
						ibd<-ibd+1
					}
				}				
				if(ibd==0){
					N0<-N0+1
				} else if (ibd==1){
					N1<-N1+1
				} else {
					N2<-N2+1
				}
				zscores<-vector()
				Ntotal<-Ntotal+1
				if (Ntotal==250){					
					p.value250<-chisq.test(c(N0,N1,N2),p=c(1,2,1),rescale.p=TRUE)$p.value
					logged250<-append(logged250,-log10(p.value250))					
				}
				if (Ntotal==500) {					
					p.value500<-chisq.test(c(N0,N1,N2),p=c(1,2,1),rescale.p=TRUE)$p.value	
					logged500<-append(logged500,-log10(p.value500))					
				}

			}			
		}
		##################################
		# After N families of simulation  #
		# I should be able to calculate   #
		# p-values from the table of IBDs #
		##################################			
		IBD<-matrix(c(N0,N1,N2),1,3)						
		p.value<-chisq.test(IBD,p=c(1,2,1),rescale.p=TRUE)$p.value	
		logged<-append(logged,-log10(p.value)) 
		r<-r+1
	} #after all replicates are done
	N<-rep(250,length(logged250))
	N<-append(N, rep(500, length(logged500)))
	N<-append(N, rep(1000,length(logged)))
	allLogs<-c(logged250,logged500,logged)	
storage<-as.matrix(cbind(N,allLogs))
write.table(storage,file=template, row.names=FALSE,col.names=FALSE)
