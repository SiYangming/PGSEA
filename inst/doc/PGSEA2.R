### R code from vignette source 'PGSEA2.Rnw'

###################################################
### code chunk number 1: PGSEA2.Rnw:36-41
###################################################
library(PGSEA)
library(GEOquery)
library(GSEABase)
gse <- getGEO("GSE7023",GSEMatrix=TRUE)
#load("gse.rda")


###################################################
### code chunk number 2: PGSEA2.Rnw:45-49
###################################################
subtype <- gsub("\\.", "_",gsub("subtype: ", "", phenoData(gse[[1]])$"characteristics_ch1"))
pheno <- new("AnnotatedDataFrame", data = data.frame(subtype), varMetadata = data.frame(labelDescription="subtype"))
rownames(pheno@data) <- colnames(exprs(gse[[1]]))
eset <- new("ExpressionSet", exprs = exprs(gse[[1]]), phenoData = pheno)


###################################################
### code chunk number 3: PGSEA2.Rnw:54-56
###################################################
data(VAIgsc)
details(VAIgsc[[1]])


###################################################
### code chunk number 4: PGSEA2.Rnw:63-64
###################################################
pg <- PGSEA(eset, VAIgsc, ref=which(subtype=="NO"))


###################################################
### code chunk number 5: PGSEA2.Rnw:69-70
###################################################
pg[5:8, 5:8]


###################################################
### code chunk number 6: PGSEA2.Rnw:78-80
###################################################
range(pg, finite=TRUE)
smcPlot(pg, col=.rwb, scale=c(-15, 15))


###################################################
### code chunk number 7: PGSEA2.Rnw:87-88
###################################################
smcPlot(pg, factor(subtype), col=.rwb, scale=c(-15, 15), margins=c(1, 1, 6, 9), show.grid=TRUE, r.cex=.75)


###################################################
### code chunk number 8: PGSEA2.Rnw:96-108
###################################################
pgNF <- PGSEA(eset, VAIgsc, ref=which(subtype=="NO"), p.value=NA)

library(limma)

design <- model.matrix(~ -1+factor(subtype))
colnames(design) <- names(table(subtype))
fit <- lmFit(pgNF, design)
contrast.matrix <- makeContrasts(P2B-NO , levels=design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

topTable(fit, n=10)[,c("logFC","t","adj.P.Val")]


###################################################
### code chunk number 9: PGSEA2.Rnw:114-115
###################################################
smcPlot(pg[as.numeric(rownames(topTable(fit, n=10))),], factor(subtype,levels=c("P1","P2B")), col=.rwb, scale=c(-15, 15), margins=c(1, 1, 6, 19), show.grid=TRUE, r.cex=.75)


###################################################
### code chunk number 10: PGSEA2.Rnw:127-137
###################################################
gos <- go2smc()
pg <- PGSEA(eset, gos, ref=which(subtype=="NO"))
pgNF <- PGSEA(eset, gos, ref=which(subtype=="NO"), p.value=NA)
#load("ABC.rda")
design <- model.matrix(~ -1+factor(subtype))
colnames(design) <- names(table(subtype))
fit <- lmFit(pgNF, design)
contrast.matrix <- makeContrasts(P2B-P1,levels=design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)


###################################################
### code chunk number 11: PGSEA2.Rnw:143-144
###################################################
smcPlot(pg[as.numeric(rownames(topTable(fit,n=30,resort.by="logFC"))),], factor(subtype,levels=c("P1","P2B")), col=.rwb, scale=c(-15, 15), margins=c(1, 1, 6, 19), show.grid=TRUE, r.cex=.75)


