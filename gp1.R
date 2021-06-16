##--------------------------------------------------------------------------------------
## Step 1: Calculate weighted average mlevel within group
## Step 2: Take the difference between groups
## Step 3: Local quadratic fit the difference, get beta0, beta1 and CIs
## Step 4: Determine region using the CI of beta1
## Step 5: Test the average beta0 within the region equals 0 or not
##-------------------------------------------------------------------------------------

##-------------------------------------------------------------------------------------
## Number of sites sequenced differs from person to person, get the shared sites first 
rm(list=ls())
datpath  <- "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
filelist <- list.files(path=datpath,pattern=".summary.rda")
load(paste(datpath,filelist[1],sep=""))
load(paste(datpath,"GSM1181642_Rep1.Proliferating.chr21.summary.rda",sep=""))
rep1 <- chr21
rm(chr21)
load(paste(datpath,"GSM1181646_Rep2.Proliferating.chr21.summary.rda",sep=""))
rep2 <- chr21
rm(chr21)
load(paste(datpath,"GSM1181647_Rep3.Proliferating.chr21.summary.rda",sep=""))
rep3 <- chr21
rm(chr21)

temp1 	<- c(rep1$pos,rep2$pos)
pos12 	<- temp1[duplicated(temp1)]
temp2 	<- c(pos12,rep3$pos)
shared.site	<- temp2[duplicated(temp2)]
save(shared.site,file=paste(datpath,"/chr21.shared.site.rda",sep=""))

##-------------------------------------------------------------------------------------
## Combine all data from six samples into one file 
## Calculate the weighted average within group and diff between groups

rm(list=ls())
datpath  <- "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
filelist <- list.files(path=datpath,pattern="chr21.summary.rda")
load(paste(datpath,filelist[1],sep=""))
load(paste(datpath,"chr21.shared.site.rda",sep="")) 
chr21_all <- chr21[,c(1,2,3,4)][which(chr21$pos %in% shared.site),]
rm(chr21)

for(ifile in 2:6){
	load(paste(datpath,filelist[ifile],sep=""))
	fil.chr21 <- chr21[which(chr21$pos %in% shared.site),]
	chr21_all <- cbind.data.frame(chr21_all,fil.chr21[,c(3,4)])
	rm(fil.chr21)
}
colnames(chr21_all) <- c("chr","pos",
						 "p1_mcount","p1_ucount",
						 "p2_mcount","p2_ucount",
						 "p3_mcount","p3_ucount",
						 "s1_mcount","s1_ucount",
						 "s2_mcount","s2_ucount",
						 "s3_mcount","s3_ucount")

chr21_all$weighted_p <- apply(chr21_all[,c(3,5,7)],1,sum)/apply(chr21_all[,3:8],1,sum) 
chr21_all$weighted_s <- apply(chr21_all[,c(9,11,13)],1,sum)/apply(chr21_all[,9:14],1,sum) 
chr21_all$weighted_diff <- chr21_all$weighted_p-chr21_all$weighted_s 
chr21_all <- chr21_all[-which(is.na(chr21_all$weighted_diff)),] # remove the NaN value

save(chr21_all,file=paste(datpath,"/chr21.all.combine.weighted.rda",sep=""))

##-------------------------------------------------------------------------------------
## Before smooth, we need to determine the widow size for each site
## Following the BSmooth idea, the widow contains at least 70 sites and at least 2kb
## For the normal kernel, the bandwidth (h) is set to be window_size/4

rm(list=ls())
library(parallel)
datpath = "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
source("/net/mulan/jiaqiang/course/biostat885/group_project/code/DMRsmooth2.R")
load(paste(datpath,"/chr21.all.combine.weighted.rda",sep=""))

## set the first site as the origin and transform into kb unit
loc 	<- (chr21_all$pos-chr21_all$pos[1])/1000 
h_chr21 <- unlist(mclapply(1:length(loc),function(x){WinCal(x,h0=0.5,loc=loc)},mc.cores=20))
save(h_chr21,file=paste(datpath,"/chr21.bandwidth.rda",sep=""))


##-------------------------------------------------------------------------------------
## Split the data into chunks for parallel computing
## 200 sites overlap at boundary 
rm(list=ls())
datpath = "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
source("/net/mulan/jiaqiang/course/biostat885/group_project/code/DMRsmooth2.R")
load(paste(datpath,"/chr21.all.combine.weighted.rda",sep=""))
load(paste(datpath,"/chr21.bandwidth.rda",sep=""))
loc <- (chr21_all$pos-chr21_all$pos[1])/1000 
## difference in percentage
wml <- chr21_all$weighted_diff*100 
fitdata <- cbind.data.frame(loc=loc,diff=wml)

chunkdata <- chunk.h <- list()
for(ichunk in 1:362){
	if(ichunk==362){
			chunkdata[[ichunk]] <- fitdata[(1+2000*(ichunk-1)):nrow(fitdata),]
			chunk.h[[ichunk]] 	<- h_chr21[(1+2000*(ichunk-1)):length(h_chr21)]
		}else{
			chunkdata[[ichunk]] <- fitdata[(1+2000*(ichunk-1)):(2000*ichunk+200),] 
			chunk.h[[ichunk]] 	<- h_chr21[(1+2000*(ichunk-1)):(2000*ichunk+200)] 
		}
}

save(chunkdata,file=paste(datpath,"/chr21.chunkdata.rda",sep=""))
save(chunk.h,file=paste(datpath,"/chr21.chunk.h.rda",sep=""))

##-------------------------------------------------------------------------------------
## LQfit 
rm(list=ls())
library(parallel)
datpath = "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
source("/net/mulan/jiaqiang/course/biostat885/group_project/code/DMRsmooth2.R")
load(paste(datpath,"/chr21.chunkdata.rda",sep=""))
load(paste(datpath,"/chr21.chunk.h.rda",sep=""))


para.func <- function(fitdata,h){
	beta    <- LQfit(fitdata, h)
	CIp   	<- AddCI(fitdata, h, beta)
	out 	<- cbind(fitdata[ , 1], round(beta, 5), CIp)
	return(out)
}

# res_chunk100 <- mclapply(1:100,function(x){para.func(chunkdata[[x]],chunk.h[[x]])},mc.cores=15)
res_chunk200 <- mclapply(101:200,function(x){para.func(chunkdata[[x]],chunk.h[[x]])},mc.cores=15)
# save(res_chunk100,file=paste(datpath,"/chr21.chunk100.res.rda",sep=""))
save(res_chunk200,file=paste(datpath,"/chr21.chunk200.res.rda",sep=""))


res_chunk300 <- mclapply(201:300,function(x){para.func(chunkdata[[x]],chunk.h[[x]])},mc.cores=15)
# save(res_chunk100,file=paste(datpath,"/chr21.chunk100.res.rda",sep=""))
save(res_chunk300,file=paste(datpath,"/chr21.chunk300.res.rda",sep=""))

res_chunk362 <- mclapply(301:362,function(x){para.func(chunkdata[[x]],chunk.h[[x]])},mc.cores=20)
# save(res_chunk100,file=paste(datpath,"/chr21.chunk100.res.rda",sep=""))
save(res_chunk362,file=paste(datpath,"/chr21.chunk362.res.rda",sep=""))



##-------------------------------------------------------------------------------------
## Combine Result
rm(list=ls())
datpath = "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
outpath = "/net/mulan/jiaqiang/course/biostat885/group_project/output/LQFIT/"
load(paste(datpath,"/chr21.chunk100.res.rda",sep=""))
load(paste(datpath,"/chr21.chunk200.res.rda",sep=""))
load(paste(datpath,"/chr21.chunk300.res.rda",sep=""))
load(paste(datpath,"/chr21.chunk362.res.rda",sep=""))

## Remove the first 100 and last 100 overlapped sites in each chunk
res_chunk100_clean <- lapply(2:100,function(x){return(res_chunk100[[x]][-c(1:100,2101:2200),])})
res_chunk200_clean <- lapply(1:100,function(x){return(res_chunk200[[x]][-c(1:100,2101:2200),])})
res_chunk300_clean <- lapply(1:100,function(x){return(res_chunk300[[x]][-c(1:100,2101:2200),])})
res_chunk362_clean <- lapply(1:61,function(x){return(res_chunk362[[x]][-c(1:100,2101:2200),])})

part1 <- do.call(rbind,res_chunk100_clean)
part2 <- do.call(rbind,res_chunk200_clean)
part3 <- do.call(rbind,res_chunk300_clean)
part4 <- do.call(rbind,res_chunk362_clean)

chr21.lqf.v1 <- rbind(res_chunk100[[1]][-c(2101:2200),],part1,part2,part3,part4,res_chunk362[[62]][-c(1:100),])
chr21.lqf.v1 <- as.data.frame(chr21.lqf.v1)
colnames(chr21.lqf.v1) <- c("loc","beta0","beta1","upperCI","lowerCI","Pvalue")
save(chr21.lqf.v1,file=paste(outpath,"/chr21.lqf.coeff.res.v1.rda",sep=""))


##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
## 400 sites overlap at boundary 
## Split the data into chunks for parallel computing 
rm(list=ls())
datpath = "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
source("/net/mulan/jiaqiang/course/biostat885/group_project/code/DMRsmooth2.R")
load(paste(datpath,"/chr21.all.combine.weighted.rda",sep=""))
load(paste(datpath,"/chr21.bandwidth.rda",sep=""))
loc <- (chr21_all$pos-chr21_all$pos[1])/1000 
## difference in percentage
wml <- chr21_all$weighted_diff*100 
fitdata <- cbind.data.frame(loc=loc,diff=wml)

temp <- fitdata[1:10200,]

chunkdata2 <- chunk.bw2 <- list()
for(ichunk in 1:362){
	if(ichunk==362){
			chunkdata2[[ichunk]] 	<- fitdata[(1+2000*(ichunk-1)):nrow(fitdata),]
			chunk.bw2[[ichunk]] 	<- h_chr21[(1+2000*(ichunk-1)):length(h_chr21)]
		}else{
			chunkdata2[[ichunk]] 	<- fitdata[(1+2000*(ichunk-1)):(2000*ichunk+400),] 
			chunk.bw2[[ichunk]] 	<- h_chr21[(1+2000*(ichunk-1)):(2000*ichunk+400)] 
		}
}

save(chunkdata2,file=paste(datpath,"/chr21.chunkdata2.rda",sep=""))
save(chunk.bw2,file=paste(datpath,"/chr21.chunk.bw2.rda",sep=""))

##------------------------------------------------------------------------------------
rm(list=ls())
library(parallel)
datpath = "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
source("/net/mulan/jiaqiang/course/biostat885/group_project/code/DMRsmooth2.R")
load(paste(datpath,"/chr21.chunkdata2.rda",sep=""))
load(paste(datpath,"/chr21.chunk.bw2.rda",sep=""))

para.func <- function(fitdata,h){
	beta    <- LQfit(fitdata, h)
	CIp   	<- AddCI(fitdata, h, beta)
	out 	<- cbind(fitdata[ , 1], round(beta, 5), CIp)
	return(out)
}

res_c100_v2 <- mclapply(1:100,function(x){para.func(chunkdata2[[x]],chunk.bw2[[x]])},mc.cores=20)
save(res_c100_v2,file=paste(datpath,"/chr21.chunk100.res.v2.rda",sep=""))
res_c200_v2 <- mclapply(101:200,function(x){para.func(chunkdata2[[x]],chunk.bw2[[x]])},mc.cores=20)
save(res_c200_v2,file=paste(datpath,"/chr21.chunk200.res.v2.rda",sep=""))
res_c300_v2 <- mclapply(201:300,function(x){para.func(chunkdata2[[x]],chunk.bw2[[x]])},mc.cores=20)
save(res_c300_v2,file=paste(datpath,"/chr21.chunk300.res.v2.rda",sep=""))
res_c362_v2 <- mclapply(301:362,function(x){para.func(chunkdata2[[x]],chunk.bw2[[x]])},mc.cores=20)
save(res_c362_v2,file=paste(datpath,"/chr21.chunk362.res.v2.rda",sep=""))

##-------------------------------------------------------------------------------------
## Combine Result
rm(list=ls())
datpath = "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
outpath = "/net/mulan/jiaqiang/course/biostat885/group_project/output/LQFIT/"
load(paste(datpath,"/chr21.chunk100.res.v2.rda",sep=""))
load(paste(datpath,"/chr21.chunk200.res.v2.rda",sep=""))
load(paste(datpath,"/chr21.chunk300.res.v2.rda",sep=""))
load(paste(datpath,"/chr21.chunk362.res.v2.rda",sep=""))

res_chunk100_clean <- lapply(2:100,function(x){return(res_c100_v2[[x]][-c(1:200,2201:2400),])})
res_chunk200_clean <- lapply(1:100,function(x){return(res_c200_v2[[x]][-c(1:200,2201:2400),])})
res_chunk300_clean <- lapply(1:100,function(x){return(res_c300_v2[[x]][-c(1:200,2201:2400),])})
res_chunk362_clean <- lapply(1:61,function(x){return(res_c362_v2[[x]][-c(1:200,2201:2400),])})

part1 <- do.call(rbind,res_chunk100_clean)
part2 <- do.call(rbind,res_chunk200_clean)
part3 <- do.call(rbind,res_chunk300_clean)
part4 <- do.call(rbind,res_chunk362_clean)

chr21.lqf.v2 <- rbind(res_c100_v2[[1]][-c(2201:2400),],part1,part2,part3,part4,res_c362_v2[[62]][-c(1:200),])
save(chr21.lqf.v2,file=paste(outpath,"/chr21.lqf.coeff.res.v2.rda",sep=""))







##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
## Overlap 1000

rm(list=ls())
datpath = "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
source("/net/mulan/jiaqiang/course/biostat885/group_project/code/DMRsmooth2.R")
load(paste(datpath,"/chr21.all.combine.weighted.rda",sep=""))
load(paste(datpath,"/chr21.bandwidth.rda",sep=""))
loc <- (chr21_all$pos-chr21_all$pos[1])/1000 
## difference in percentage
wml <- chr21_all$weighted_diff*100 
fitdata <- cbind.data.frame(loc=loc,diff=wml)

chunkdata3 <- chunk.bw3 <- list()
for(ichunk in 1:723){
	if(ichunk==723){
			chunkdata3[[ichunk]] 	<- fitdata[(1+1000*(ichunk-1)):nrow(fitdata),]
			chunk.bw3[[ichunk]] 	<- h_chr21[(1+1000*(ichunk-1)):length(h_chr21)]
		}else{
			chunkdata3[[ichunk]] 	<- fitdata[(1+1000*(ichunk-1)):(1000*ichunk+1000),] 
			chunk.bw3[[ichunk]] 	<- h_chr21[(1+1000*(ichunk-1)):(1000*ichunk+1000)] 
		}
}

save(chunkdata3,file=paste(datpath,"/chr21.chunkdata3.rda",sep=""))
save(chunk.bw3,file=paste(datpath,"/chr21.chunk.bw3.rda",sep=""))




##------------------------------------------------------------------------------------
rm(list=ls())
library(parallel)
datpath = "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
source("/net/mulan/jiaqiang/course/biostat885/group_project/code/DMRsmooth2.R")
load(paste(datpath,"/chr21.chunkdata3.rda",sep=""))
load(paste(datpath,"/chr21.chunk.bw3.rda",sep=""))

para.func <- function(fitdata,h){
	beta    <- LQfit(fitdata, h)
	CIp   	<- AddCI(fitdata, h, beta)
	out 	<- cbind(fitdata[ , 1], round(beta, 5), CIp)
	return(out)
}

res_c100_v3 <- mclapply(1:100,function(x){para.func(chunkdata3[[x]],chunk.bw3[[x]])},mc.cores=20)
save(res_c100_v3,file=paste(datpath,"/chr21.chunk100.res.v3.rda",sep=""))

res_c200_v3 <- mclapply(101:200,function(x){para.func(chunkdata3[[x]],chunk.bw3[[x]])},mc.cores=20)
save(res_c200_v3,file=paste(datpath,"/chr21.chunk200.res.v3.rda",sep=""))

res_c300_v3 <- mclapply(201:300,function(x){para.func(chunkdata3[[x]],chunk.bw3[[x]])},mc.cores=20)
save(res_c300_v3,file=paste(datpath,"/chr21.chunk300.res.v3.rda",sep=""))

res_c400_v3 <- mclapply(301:400,function(x){para.func(chunkdata3[[x]],chunk.bw3[[x]])},mc.cores=20)
save(res_c400_v3,file=paste(datpath,"/chr21.chunk400.res.v3.rda",sep=""))

res_c500_v3 <- mclapply(401:500,function(x){para.func(chunkdata3[[x]],chunk.bw3[[x]])},mc.cores=20)
save(res_c500_v3,file=paste(datpath,"/chr21.chunk500.res.v3.rda",sep=""))

res_c600_v3 <- mclapply(501:600,function(x){para.func(chunkdata3[[x]],chunk.bw3[[x]])},mc.cores=20)
save(res_c600_v3,file=paste(datpath,"/chr21.chunk600.res.v3.rda",sep=""))

res_c700_v3 <- mclapply(601:700,function(x){para.func(chunkdata3[[x]],chunk.bw3[[x]])},mc.cores=20)
save(res_c700_v3,file=paste(datpath,"/chr21.chunk700.res.v3.rda",sep=""))

res_c723_v3 <- mclapply(701:723,function(x){para.func(chunkdata3[[x]],chunk.bw3[[x]])},mc.cores=20)
save(res_c723_v3,file=paste(datpath,"/chr21.chunk723.res.v3.rda",sep=""))



# res_c724_v3 <- mclapply(721:723,function(x){para.func(chunkdata3[[x]],chunk.bw3[[x]])},mc.cores=20)
# save(res_c724_v3,file=paste(datpath,"/chr21.chunk724.res.v3.rda",sep=""))

# para.func(chunkdata3[[722]],chunk.bw3[[722]])

##-------------------------------------------------------------------------------------
## Combine Result
rm(list=ls())
datpath = "/net/mulan/jiaqiang/course/biostat885/group_project/data/senescent/"
outpath = "/net/mulan/jiaqiang/course/biostat885/group_project/output/LQFIT/"
load(paste(datpath,"/chr21.chunk100.res.v3.rda",sep=""))
load(paste(datpath,"/chr21.chunk200.res.v3.rda",sep=""))
load(paste(datpath,"/chr21.chunk300.res.v3.rda",sep=""))
load(paste(datpath,"/chr21.chunk400.res.v3.rda",sep=""))
load(paste(datpath,"/chr21.chunk500.res.v3.rda",sep=""))
load(paste(datpath,"/chr21.chunk600.res.v3.rda",sep=""))
load(paste(datpath,"/chr21.chunk700.res.v3.rda",sep=""))
load(paste(datpath,"/chr21.chunk723.res.v3.rda",sep=""))

res_chunk100_clean <- lapply(2:100,function(x){return(res_c100_v3[[x]][-c(1:500,1501:2000),])})
res_chunk200_clean <- lapply(1:100,function(x){return(res_c200_v3[[x]][-c(1:500,1501:2000),])})
res_chunk300_clean <- lapply(1:100,function(x){return(res_c300_v3[[x]][-c(1:500,1501:2000),])})
res_chunk400_clean <- lapply(1:100,function(x){return(res_c400_v3[[x]][-c(1:500,1501:2000),])})
res_chunk500_clean <- lapply(1:100,function(x){return(res_c500_v3[[x]][-c(1:500,1501:2000),])})
res_chunk600_clean <- lapply(1:100,function(x){return(res_c600_v3[[x]][-c(1:500,1501:2000),])})
res_chunk700_clean <- lapply(1:100,function(x){return(res_c700_v3[[x]][-c(1:500,1501:2000),])})
res_chunk723_clean <- lapply(1:22,function(x){return(res_c723_v3[[x]][-c(1:500,1501:2000),])})

part1 <- do.call(rbind,res_chunk100_clean)
part2 <- do.call(rbind,res_chunk200_clean)
part3 <- do.call(rbind,res_chunk300_clean)
part4 <- do.call(rbind,res_chunk400_clean)
part5 <- do.call(rbind,res_chunk500_clean)
part6 <- do.call(rbind,res_chunk600_clean)
part7 <- do.call(rbind,res_chunk700_clean)
part8 <- do.call(rbind,res_chunk723_clean)

chr21.lqf.v3 <- rbind(res_c100_v3[[1]][-c(1501:2000),],
						part1,part2,part3,part4,part5,part6,part7,part8,
						res_c723_v3[[23]][-c(1:500),])

save(chr21.lqf.v3,file=paste(outpath,"/chr21.lqf.coeff.res.v3.rda",sep=""))


##-----------------------------
# rm(list=ls())
# out1 <- read.table("/net/mulan/jiaqiang/Online_Data/GEO/senescent_cell/GSE48580_Proliferating.vs.Senescent.PooledG.DMR.up.bed")
# res1  <- out1[out1[,1]=="chr18",]

# out2 <- read.table("/net/mulan/jiaqiang/Online_Data/GEO/senescent_cell/GSE48580_Proliferating.vs.Senescent.PooledG.DMR.down.bed")
# res2  <- out2[out2[,1]=="chr18",]















