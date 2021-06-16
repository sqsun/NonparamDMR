load("/net/wonderland/home/shiquans/temp/course/biostat885/group_project/output/LQFIT/chr21.lqf.coeff.res.v3.rda")
library(dplyr)
library(scales)

mat <- data.frame(chr21.lqf.v3[11000:723000,])
#mat <- data.frame(chr21.lqf.v3[11000:21000,])
names(mat)[1] <- "loc"
mat2 <- mat; mat2$beta0 <- abs(mat$beta0)

###### Fine CI windows ######
IsEP <- function(dat, idx, nb) {
	re <- FALSE
	if (dat$beta0[idx]>0) {
	if (dat$beta1[idx-1]>0 & dat$beta1[idx+1]<0) {
		if (dat$upperCI[idx]>0 & dat$lowerCI[idx]<0 & (dat$beta0[idx]==max(dat$beta0[(idx-nb):(idx+nb)]))) {re <- TRUE}
	}
	else if (dat$beta1[idx-1]<0 & dat$beta1[idx]<0 & dat$beta1[idx+1]<0) {
		if (dat$upperCI[idx]>0 & dat$lowerCI[idx]<0 & (dat$beta0[idx]==max(dat$beta0[(idx-nb):(idx+nb)]))) {re <- TRUE}
	}
	else if (dat$beta1[idx-1]>0 & dat$beta1[idx]>0 & dat$beta1[idx+1]>0) {
		if (dat$upperCI[idx]>0 & dat$lowerCI[idx]<0 & (dat$beta0[idx]==max(dat$beta0[(idx-nb):(idx+nb)]))) {re <- TRUE}
	}
	}
	return(re)
}

GetEP <- function(dat, nb) {
	idx.vec <- (nb+1):(nrow(dat)-nb)
	mark.vec <- sapply(idx.vec, IsEP, dat=dat, nb=nb)
	mark.vec <- c(rep(FALSE, nb), mark.vec, rep(FALSE, nb))
	return(mark.vec)
}

GetCI <- function(dat, idx, nb) {
#print(idx)
	upperCI <- dat$upperCI[idx]; lowerCI <- dat$lowerCI[idx]
	u <- d <- 1
	while(dat$beta1[idx-u]<upperCI & u<(nb+1)) {u <- u+1}
	while(dat$beta1[idx+d]>lowerCI & d<(nb+1)) {d <- d+1}
	if (max(u,d)<=nb) {re <- c(ind=TRUE, idx=idx, u=u, d=d, loc=dat$loc[idx], beta0=dat$beta0[idx], lower=dat$loc[idx-u], upper=dat$loc[idx+d])}
	else if (min(u,d)<=nb) {u<-d<-min(u,d); re <- c(ind=TRUE, idx=idx, u=u, d=d, loc=dat$loc[idx], beta0=dat$beta0[idx], lower=dat$loc[idx-u], upper=dat$loc[idx+d])}
	else {re <- c(ind=FALSE, idx=idx, u=NA, d=NA, loc=NA, beta0=NA, lower=NA, upper=NA)}
	
	return(re)
}

GetGoodCI <- function(dat, nb) {
	peakidx <- (1:nrow(dat))[GetEP(dat, nb)]
	
	CI.result <- data.frame(t(sapply(peakidx, GetCI, dat=dat, nb=nb)))
	re <- CI.result[CI.result$ind==1,]
	return(re)
}

CI.mat2 <- GetGoodCI(dat=mat2, nb=10); dim(CI.mat2)

CI.mat <- CI.mat2
positive.idx <- mat$beta0[CI.mat2$idx]/CI.mat2$beta0
CI.mat$beta0 <- mat$beta0[CI.mat2$idx]*positive.idx

write.table(CI.mat, file="results_CI_2018-04-16.txt", col.names=TRUE)

###### Find Islands ######
FindIsland <- function(dat, CI.mat, idx, distcut, sitenum) {
	#print(idx)
	repnum <- sitenum-2
	u <- CI.mat$u[CI.mat$idx==idx]; d <- CI.mat$d[CI.mat$idx==idx]; idx0 <- idx-u
	
	if (u+d+1<sitenum) {out <- rep(NA,4)}
	else {
	submat <- mat[(idx-u):(idx+d),]
	loc2 <- lag(submat$loc)
	diffs <- (submat$loc-loc2)[-1]
	diffs2 <- which(diffs<distcut)
	
	tmp <- rle(diff(diffs2))
	s1 <- diffs2[1]
	re <- NULL
	for (i in 1:length(tmp$values)) {
		if (tmp$values[i]==1 & tmp$lengths[i]>=repnum) {startidx <- s1} else {startidx <- NA}
		s1 <- s1 + (tmp$lengths[i])*(tmp$values[i])
		#print(s1)
		if (tmp$values[i]==1 & tmp$lengths[i]>=repnum) {endidx <- s1} else {endidx <- NA}
		re <- rbind(re, c(startidx=startidx, endidx=endidx))
		#print(re)
	}
	re <- re[!is.na(re[,1]),]
	#print(re)
	#print(nrow(re))
	
	if (is.null(re)) {out <- rep(NA, 4)} 
	else if (length(re)==0) {out <- rep(NA, 4)}
	else {
	re <- matrix(re, ncol=2)
	#withinpeakid <- rep(1:nrow(re), (re[,2]-re[,1])+2, each=T)
	#peakidx <- rep(idx, length(withinpeakid))
	peakidx <- rep(idx, nrow(re))
	peakloc <- rep(dat$loc[idx], nrow(re))
	startloc <- submat$loc[re[,1]]
	endloc <- submat$loc[re[,2]+1]
	
	out <- t(cbind(peakidx, peakloc, startloc, endloc))
	}
	}
	return(out)
	
}

FindAllIslands <- function(dat, CI.mat, distcut, sitenum) {
	tmp <- sapply(CI.mat$idx, FindIsland, dat=mat, CI.mat=CI.mat, 0.3, sitenum)
	re <- na.omit(data.frame(matrix(unlist(tmp), ncol=4, byrow=T)))
	names(re) <- c("peakidx", "peakloc", "startloc", "endloc")	
	return(re)
}

island.mat <- FindAllIslands(dat=mat, CI.mat, distcut=0.3, sitenum=10)
table(island.mat$peakidx); nrow(island.mat) #2334; 1085 

Do.onettest <- function(dat, island.mat, island.mat.idx) {
	startloc <- island.mat$startloc[island.mat.idx]
	endloc <- island.mat$endloc[island.mat.idx]
	
	setdat <- dat[(match(startloc, dat$loc)):(match(endloc, dat$loc)),]
	pval <- t.test(setdat$beta0)$p.value
	average.beta0 <- mean(setdat$beta0)
	return(c(pval=pval, average.beta0=average.beta0))
}

Do.allttest <- function(dat, island.mat) {
	ttest.pval <- t(sapply(1:nrow(island.mat), Do.onettest, dat=dat, island.mat=island.mat))
	re <- cbind(island.mat, ttest.pval)
	return(re)
}

island.re <- Do.allttest(dat=mat, island.mat); dim(island.re) #1085 islands

write.table(island.re, file="results_island_2018-04-16.txt", col.names=TRUE)

table(island.re$pval<(0.05/nrow(island.re))) #1077 significant islands

###### Make plots
####Plot the results
sub.mat <- mat[23000:25500,]
x <- sub.mat$loc
y <- sub.mat$beta0
l <- length(x)
ma <- max(x)
beta0 <- sub.mat$beta0
beta1 <- sub.mat$beta1
upperCI <- sub.mat$upperCI
lowerCI <- sub.mat$lowerCI

sub.CI.mat <- CI.mat[CI.mat$loc %in% sub.mat$loc,]
sub.island.re <- island.re[(island.re$peakloc %in% sub.mat$loc) & (island.re$pval<(0.05/nrow(island.re))),]

## nonparametric fitted curve with 95% CI
par(mfrow=c(3,1))

# if plotting CI windows
plot( x, y, ylim=c(-2, max(y)+2), type = 'l', lwd=2,  xlab = 'Chromosomal coordinates (in kb)', ylab = expression(beta[0]), main="Nonparametric Fitted Curve with 95% CI Windows",cex.axis=1.6, cex.lab=1.8, cex.main=2)
vx<-as.vector( sub.CI.mat$loc)
ll <- length( vx)
vx1<-as.vector( sub.CI.mat$lower)
vx2<-as.vector( sub.CI.mat$upper)
vy <- as.vector(sub.CI.mat$beta0)

points( vx, vy+1, pch = 25, col=1, bg=1)
for ( i in 1:ll){
       x1 <- vx1[i]
       x2 <- vx2[i]
       y1 <- vy[i]
       #polygon( c( x1, x1, x2, x2), c( 0, 0.7*y1, 0.7*y1, 0), col = 16)
       polygon( c( x1, x1, x2, x2), c( 0, y1-3, y1-3, 0), col = "gray50", border=FALSE)       
}

# if plotting islands
plot( x, y, ylim=c(-5, max(y)+2), type = 'l', lwd=2,  xlab = 'Chromosomal coordinates (in kb)', ylab = expression(beta[0]), main="Nonparametric Fitted Curve with Differentially Expressed DMR Islands",cex.axis=1.6, cex.lab=1.8, cex.main=2)
vx<-as.vector( sub.island.re$peakloc)
ll <- length( vx)
vx1<-as.vector( sub.island.re$startloc)
vx2<-as.vector( sub.island.re$endloc)
vy <- as.vector(sub.island.re$average.beta0)

points( vx, vy+2, pch = 25, col="darkslateblue", bg=1)
for ( i in 1:ll){
       x1 <- vx1[i]
       x2 <- vx2[i]
       y1 <- vy[i]
       #polygon( c( x1, x1, x2, x2), c( 0, 0.7*y1, 0.7*y1, 0), col = 16)
       polygon( c( x1, x1, x2, x2), c( 0, y1-3, y1-3, 0), col=alpha("darkslateblue", 0.5), border=FALSE)
}


##first-derivative curve with 95% CI
plot( x, beta1, ylim=c(min(lowerCI)-2, max(upperCI+2)), type = 'l', xlab = 'Chromosomal coordinates (in kb)', main='Estimated First Derivative Curve with 95% Confidence Interval', ylab = expression(beta[1]),cex.axis=1.6, cex.lab=1.8, cex.main=2)
polygon(c(x, rev(x)), c(lowerCI, rev(upperCI)), col=alpha("darkorange", 0.5), border=FALSE)
abline(h=0)
points(vx, rep(0, length(vx)), pch = "x")

