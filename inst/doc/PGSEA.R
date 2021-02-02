### R code from vignette source 'PGSEA.Rnw'

###################################################
### code chunk number 1: PGSEA.Rnw:40-43
###################################################
	library(PGSEA)
	basic <- new("smc",ids=c("gene a","gene b"),reference="simple smc")
	str(basic)


###################################################
### code chunk number 2: PGSEA.Rnw:50-53
###################################################
	datadir <- system.file("extdata", package = "PGSEA")
	sample <- readGmt(file.path(datadir, "sample.gmt"))
	str(sample[1])


###################################################
### code chunk number 3: PGSEA.Rnw:58-61
###################################################

	data(nbEset)
	pg <- PGSEA(nbEset,cl=sample,ref=1:5)


###################################################
### code chunk number 4: PGSEA.Rnw:66-69
###################################################
	sub <- factor(c(rep(NA,5),rep("NeuroB",5),rep("NeuroB_MYC+",5)))
	smcPlot(pg[,],sub,scale=c(-12,12),show.grid=T,margins=c(1,1,7,13),col=.rwb)
	


###################################################
### code chunk number 5: PGSEA.Rnw:77-80
###################################################
	mcs <- go2smc()[1:10]
	pg <- PGSEA(nbEset,cl=mcs,ref=1:5)
	


###################################################
### code chunk number 6: PGSEA.Rnw:85-88
###################################################

	smcPlot(pg[,],sub,scale=c(-12,12),show.grid=T,margins=c(1,1,7,20),col=.rwb)
	


###################################################
### code chunk number 7: PGSEA.Rnw:98-102
###################################################
	#data(VAImcs)
	data(VAIgsc)
  pg <- PGSEA(nbEset,cl=VAIgsc,ref=1:5)
	


###################################################
### code chunk number 8: PGSEA.Rnw:107-110
###################################################

	smcPlot(pg[,],sub,scale=c(-5,5),show.grid=T,margins=c(1,1,8,14),col=.rwb,r.cex=.7)
	


