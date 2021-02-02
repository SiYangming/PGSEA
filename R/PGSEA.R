
setClass("smc",representation(reference="character",desc="character",source="character",design="character",identifier="character",species="character",data="character",private="character", creator="character",ids="vector"))
  
 
  
  
PGSEA <- function (exprs, cl, range=c(25,500), ref=NULL,center=TRUE, p.value=0.005, weighted=TRUE, enforceRange=TRUE,...) {

  if (is(exprs, "ExpressionSet"))
    exprs <- exprs(exprs)

  if(!is.list(cl))
    stop("cl need to be a list")  
  if(!is.null(ref)) {
    if(!is.numeric(ref))
      stop("column index's required")
  }
  if(!is.null(ref)){
    if(options()$verbose)
      cat("Creating ratios...","\n")
    ref_mean <- apply(exprs[,ref],1,mean,na.rm=TRUE)
    exprs <- sweep(exprs,1, ref_mean,"-")
  }
  if(center)
    exprs <- scale(exprs,scale=FALSE)
  results <- matrix(NA,length(cl),ncol(exprs))
  rownames(results) <- names(cl)
  colnames(results) <- colnames(exprs)
  mode(results) <- "numeric"
  if(is.logical(p.value))
    p.results <- results

  for (i in 1:length(cl)) {
    if(class(cl[[i]]) == "smc") {
      clids <- cl[[i]]@ids
    } else if(class(cl[[i]]) %in% c("GeneColorSet","GeneSet")) {
      clids <- cl[[i]]@geneIds
    } else  {
      clids <- cl[[i]]
    }
    if(options()$verbose)
      cat("Testing region ", i, "\n")
    ix <- match(clids,rownames(exprs))
    ix <- unique(ix[!is.na(ix)])
    present <- sum(!is.na(ix))
    if(present < range[1]) {
      if(options()$verbose) cat("Skipping region ",i," because too small-",present,",\n")
      next
    }
    if(present > range[2]) {
       if(options()$verbose) cat("Skipping region ",i," because too large-",present,"\n")
      next
    }

    texprs <- exprs[ix,]
    if(any(is.na(texprs))) cat("Warning - 'NA' values within expression data, enrichment scores are estimates only.\n")
    if(!is.matrix(texprs)) texprs <-as.matrix(texprs)

	if(!weighted) stat <- try(apply(texprs, 2, t.test,...))
	else {
		try(mod <- (length(ix) ^ (1/2)) / apply(exprs, 2, sd, na.rm=TRUE))
		try(m <- apply(texprs, 2, mean,na.rm=TRUE) - apply(exprs,2,mean,na.rm=TRUE))
		stat2 <- m * mod
    p.val <- 2*pnorm(-abs(stat2))
		stat <- list()
		for(q in 1:length(stat2)){
			stat[[q]] <- list(statistic=stat2[q],p.value=p.val[q])
		}
		names(stat) <- names(stat2)
		
	}
    if (is.list(stat)) {
      ps <- unlist(lapply(stat, function(x) x$p.value))
      stat <- unlist(lapply(stat, function(x) x$statistic))
      if (!is.na(p.value)) {
        if(is.numeric(p.value)) {
          stat[ps > p.value] <- NA
        } else {
          p.results[i,] <- ps
        }
      }
    }
    results[i,] <- as.numeric(stat)
   if(enforceRange) {
    for(w in 1:ncol(texprs)) {
      if(sum(!is.na(texprs[,w])) < range[1] | sum(!is.na(texprs[,w])) > range[2] ) 
      results[i,w] <- NA
    }
   }
  }
 if(is.logical(p.value) & !is.na(p.value)) return(list(results=results,p.results=p.results))
 return(results)
}




aggregateExprs <- function(x,package="hgu133plus2",using="ENTREZID",FUN,...) {
  if(class(x) != "matrix" && !is(x, "ExpressionSet"))
    stop("need matrix or ExpressionSet")
  if(is.null(package) && is(x, "ExpressionSet"))
    package <- annotation(x)
  if(is.null(package))
    stop("annotation package name is required")
  if(!require(package,character.only=TRUE))
    stop(package," is not available")

  pPos <- paste("package",package,sep=":")
  if(length(grep(".*db",package))) package <- gsub(".db","",package)
  
  nEnv <- paste(package,using,sep="")
  Env <- get(nEnv,pos=pPos)
  if(is(x, "ExpressionSet")) {
    ids <- featureNames(x)
    x <- exprs(x)
  } else {
    ids <- rownames(x)
  }
  lls <- mget(ids,env=Env,ifnotfound=NA)
  if(length(lls)!=length(unlist(lls))) for(i in 1:length(lls)) lls[[i]] <- lls[[i]][1]
  lls <- unlist(lls)
  f <- factor(lls)
  undupx <- aggregate(x,by=list(f),FUN,...)
  rownames(undupx) <- as.character(undupx[,1])
  undupx <- as.matrix(undupx[,-1])
  return(undupx)
}

readGmt <- function(fname) {
  f <- readLines(fname)
  mc <- list()
  for(i in 1:length(f)) {
    dat <- unlist(strsplit(f[i],"\t",fixed=TRUE))
    m <- new("smc")
    m@reference <- dat[1]
    if(dat[2] != "NA") m@desc <- dat[2]
	else m@desc <- ""
    ids <- dat[3:length(dat)]
    m@ids <- ids[!(ids=="NA")]
    mc <- c(mc,list(m))
  }
  names(mc) <- unlist(lapply(mc,function(x) paste(x@reference,x@desc)))
  return(mc)
}

writeGmt <- function(fname,cl) {

  sink(fname)
  for(i in 1:length(cl)){
	cat(cl[[i]]@reference,"\t")
	cat(cl[[i]]@desc,"\t")
	cat(cl[[i]]@ids,sep="\t")
	cat("\n")
  }
  sink()
	
}


readSmc <- function(files){

  smc.list <- list()
  smc.names <- c()
  for(i in 1:length(files)){
    f <- readLines(files[i])
    header <- unlist(strsplit(f[1],"|",fixed=TRUE))
    ref <- sub(">","",header[1])
    ref <- gsub(" ","",ref)
    table.num <- try(as.numeric(header[2]))
    chip <- gsub(" ","",header[3])
    creator <- gsub(" ","",header[4])
    desc <- header[5]
    genes <- unlist(strsplit(f[2],",",fixed=TRUE))
    smc.list <- c(smc.list,list(new("smc",reference=ref,table=table.num,chip=chip,creator=creator,desc=desc,ids=genes,p.value=vector(),statistic=vector())))
    smc.names <- c(smc.names,paste(ref,table.num,sep="."))
  }
  names(smc.list) <- smc.names
  return(smc.list)
}





writeSmc <-function(x){
  if(is.null(x@desc) | is.null(x@chip) | is.null(x@creator))
    stop("A description, chip, and creator is required")
  if(length(x@reference)==0 ) {
    warning("No reference: creating 'random' unique id")
    x@reference <- paste(sample(c(0:9,letters,LETTERS),40,replace=TRUE),collapse="")
  } else {
    #stop("A references is required. Cannot autogenerate")
  }
  if(length(x@table) ==0) x@table <- 0  
  fileName <- paste(x@reference,"-",x@table,".txt",sep="")
  sink(fileName)
  cat(">",x@reference,"|",x@table,"|",x@chip,"|",x@creator, "|", x@desc,"\n")
  write.table(x@ids,sep=",",row.names=FALSE,eol=",",col.names=F,quote=FALSE)
  sink()
  return(fileName)
}


editSmc <- function (smcList,attName="creator",newAtt="changed!!"){
	if(length(newAtt)==1)	for(i in 1:length(smcList)) attr(smcList[[i]],attName) <- newAtt
	else for(i in 1:length(smcList)) attr(smcList[[i]],attName) <- newAtt[i]

	return(smcList)
	
}


scanSmc <- function(smcList,scanSlot="private",scanFor="no"){
	ix <- vector()
	for(i in 1:length(smcList)){
		if(length(grep(scanFor,attr(smcList[[i]],scanSlot)) == 1)){
			ix <- c(ix,i)
		}
	}
	return(smcList[ix])
}


convertSmc <- function(mcs,fromSpecies="h", toSpecies="r",hgX="./homologene.data"){

	fromSpecies <- switch(fromSpecies, h =9606, r = 10116, m = 10090,c=9615)
	toSpecies <- switch(toSpecies, h =9606, r = 10116, m = 10090,c=9615)

	hgX <- read.delim(hgX,col.names=c("HID","TaxID","GeneID","Symbol" ,"Protein","Protein-accession",sep="|"))[,1:3]
	
	hs <- which(hgX$TaxID %in% fromSpecies)
	rn <- which(hgX$TaxID %in% toSpecies)
	from <- hgX[hs,]
	to <- hgX[rn,]
	
	for(z in 1:length(mcs)){
		matched <- mcs[[z]]@ids
		ids <- to[match(from[match(matched,from[,"GeneID"]),"HID"],to[,"HID"]),"GeneID"]
		mcs[[z]]@ids <- ids[!is.na(ids)]
		if(z %% 50 == 0) cat("finished ",z,"of",length(mcs),"\n")
	}
	
	return(mcs)
	
}

smcPlot <- function(m,ff=NULL,skip="NO",scale=c(-3,3),na.color=par("bg"),margins=NULL,r.cex=NULL,c.cex=NULL,show.grid=F,cnames=TRUE,rnames=TRUE,grid.lty=3,clust=FALSE,...){
  hold <- as.matrix(m[nrow(m):1,])
  colnames(hold) <- colnames(m)
  m <- hold
  rm(hold)
  di <- dim(m)
  nr <- di[1]
  nc <- di[2]


if(!is.null(ff)){
  new.ix <- rep(NA,length(!is.na(ff)))
  start <- 1
  end <- NA
  for(i in levels(ff)) {
    if (i %in% skip) next;
    ix <- which(ff==i)
    if(length(ix) == 0) warning(i," is not found")
    end <- start + length(ix) -1
    new.ix[start:end] <- ix
    start <- end + 1
  }
  new.ix <- new.ix[!is.na(new.ix)]
  m <- m[,new.ix]
  ff <- ff[new.ix]

  cn <- colnames(m)
  for(i in levels(ff)) {
    ix <- which(ff==i)
    if(length(ix) == 0) next
    colnames <- rep(NA,length(ix))
    colnames[ceiling(length(ix)/2)] <- i
    cn[ix] <- colnames
  }
  colnames(m) <- cn
}

  #want to cluster?
if(clust){
    #also order based on a factor?
  if(!is.null(ff)){
    cn <- colnames(m)
    for(i in levels(ff)) {
      ix <- which(ff==i)
      if(length(ix) > 3){  
        clust.gx <- m[,ix]
        hc <- hclust(dist(t(clust.gx),method="euc"))
        m[,ix] <- m[,ix[hc$order]]
      }
    
    }      
    colnames(m) <- cn
  }
    #no factor ordering
  if(is.null(ff)){
      hc <- hclust(dist(t(m),method="euc"))
      m[,] <- m[,hc$order]
  }
}

  if(!is.null(scale)) {
    if(length(scale) == 2) {
      m[m < scale[1]] <- scale[1]
      m[m > scale[2]] <- scale[2]      
      GS <- rep(NA,length=ncol(m))
      gs <- scale[1]:scale[2]
      #GS[1:length(gs)] <- gs
      GS[1] <- scale[1]
      GS[2] <- scale[2]
    } else {
      m[m < min(scale)] <- min(scale)
      m[m > max(scale)] <- max(scale)
      GS <- rep(0,length=ncol(m))
      GS[1:length(scale)] <- scale
    }
    m <- rbind(GS,m)
  }
  
  op <- par(no.readonly = TRUE)
  if(is.null(r.cex))
    r.cex <- op$cex.axis
  if(is.null(c.cex))
    c.cex <- op$cex.axis
  if(is.null(margins)) {
    par(mar = c(1,1,3,10))
  } else {
    par(mar=margins)
   }
  image(1:ncol(m), 1:nrow(m), t(m), axes = FALSE, xlim = c(0.5,ncol(m) + 0.5), ylim = c(0.5, nrow(m) + 0.5), xlab = ",",ylab = ",",...)
  if (!is.na(na.color) & any(is.na(m))) {
    na.m <- ifelse(is.na(m), 1, NA)
    image(1:ncol(m), 1:nrow(m), t(na.m), axes = FALSE, , 
          xlab = ",", ylab = ",", col = na.color, add = T)
  }
  if(show.grid) grid(ncol(m),nrow(m),col="slategrey",lty=grid.lty)
    
  if(cnames==TRUE) axis(3, 1:ncol(m), las = 2, line = -0.5, tick = 0, labels = colnames(m), cex.axis = c.cex)
    else if(is.character(cnames)) axis(3, 1:ncol(m), las = 2, line = -0.5, tick = 0, labels = cnames, cex.axis = c.cex)
  if(rnames==TRUE) axis(4, 1:nrow(m), las = 2, line = -0.5, tick = 0, labels = rownames(m), cex.axis = r.cex)
    else if(is.character(rnames)) axis(4, 1:nrow(m), las = 2, line = -0.5, tick = 0, labels = rnames, cex.axis = r.cex)



  if(!is.null(ff)) {
	locations <- cumsum(table(ff))[1:(length(levels(ff))-1)] + 0.5
	abline(v=locations)
  }
  box()
  par(op)
}


kegg2smc <-function (min = 1, max = 284,organism="human"){
    if(organism!="human") stop("Only human for now..")
    if (!require(KEGG.db))
        stop("library(KEGG.db) is required...")
    else {
        library(KEGG.db)
        terms <- mget(ls(KEGGPATHNAME2ID), KEGGPATHNAME2ID, ifnotfound = NA)
            #necessary to add "hsa" prefix back to KEGG IDs otherwise mget fails to match character string
      	splice.id <- lapply(terms,function(x) {paste("hsa", x, sep="")})
        lists <- mget(as.character(splice.id), KEGGPATHID2EXTID, ifnotfound = NA)
        lengths <- unlist(lapply(lists, length))
        ix <- which(lengths > min & lengths < max)
        terms <- terms[ix]
        lists <- lists[ix]
        keggList <- list()
      	keggNames <- names(terms)

        for (i in 1:length(keggNames)) keggList[[i]] <- new("smc",
            reference = terms[[i]], desc = names(terms[i]),
            source = "KEGG", identifier = terms[[i]], species = "human",
            design = "ONTO", ids = lists[[i]])
        names(keggList) <- unlist(lapply(keggList, function(x) paste(x@reference,
            x@desc)))
        return(keggList)
    }
}


go2smc <- function(min=50,max=200,organism="human"){
  if(organism!="human") stop("Only human for now..")
	if(!require(org.Hs.eg.db)) stop("library(org.Hs.eg.db) is required...")
	if(!require(GO.db)) stop("library(GO.db) is required...")
	else {
		library(GO.db)
		library(org.Hs.eg.db)
#		terms <- mget(ls(GOTERM),GOTERM@datacache,ifnotfound=NA)
#		lists <- mget(names(terms),org.Hs.egGO2ALLEGS@datacache,ifnotfound=NA)
		terms <- mget(ls(GOTERM),GOTERM,ifnotfound=NA)
  	lists <- mget(names(terms),org.Hs.egGO2ALLEGS,ifnotfound=NA)
		lengths <- unlist(lapply(lists,length))
		
		ix <- which(lengths > min & lengths < max)
		
		terms <- terms[ix]
		lists <- lists[ix]
		goList <- list()
		for(i in 1:length(terms)) 
			goList[[i]] <- new("smc",
				reference=terms[[i]]@GOID,
				desc=terms[[i]]@Term,
				source="GO",
				identifier=terms[[i]]@GOID,
				species="human",
				design="ONTO",
				ids=as.character(lists[[i]])) 
		names(goList) <- unlist(lapply(goList,function(x) paste(x@reference,x@desc)))
		return(goList)
	}
}

.rwb <- c("#0000FF","#0B0BFF","#1515FF","#2020FF","#2B2BFF","#3535FF","#4040FF"
,"#4A4AFF","#5555FF","#6060FF","#6A6AFF","#7575FF","#8080FF","#8A8AFF"
,"#9595FF","#9F9FFF","#AAAAFF","#B5B5FF","#BFBFFF","#CACAFF","#D4D4FF"
,"#DFDFFF","#EAEAFF","#F4F4FF","#FFFFFF","#FFFFFF","#FFF4F4","#FFEAEA"
,"#FFDFDF","#FFD5D5","#FFCACA","#FFBFBF","#FFB5B5","#FFAAAA","#FF9F9F"
,"#FF9595","#FF8A8A","#FF8080","#FF7575","#FF6A6A","#FF6060","#FF5555"
,"#FF4A4A","#FF4040","#FF3535","#FF2B2B","#FF2020","#FF1515","#FF0B0B"
,"#FF0000")

