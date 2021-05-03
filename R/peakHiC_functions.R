extractReads <- function( fragID, vpS, vpZoom=c(1e6,1e6), reads, frags, k=31, vpSize=31 ){
  
  fragsChr <- frags
  # first we search for the fragment containing the viewpoint
  vpIdx <- match(fragID,fragsChr$fragID)
  
  # if the fragID doesn't exist, return NULL
  if(is.na(vpIdx)) {
    
    return(list( GR=NULL, vpGR=NULL, vpCov=0, plotCov=0 ))
    
  }
  
  # then we select all the fragments around the viewopint in the range provided
  vpF <- fragsChr[ queryHits( findOverlaps( ranges( fragsChr ), IRanges( start=vpS-vpZoom[1], end=vpS+vpZoom[2] ) ) ) ]
  
  if(length(vpF) < 5*vpSize){
    
    return(list( GR=NULL, vpGR=NULL, vpCov=0, plotCov=0 ))
  }
  
  # fragments pooled as viewpoint
  metaFrags <- (vpIdx-floor(vpSize/2)):(vpIdx+floor(vpSize/2))
  
  # check for beginning of chr
  if(vpIdx < vpSize) {metaFrags <- 1:vpSize}
  metaFrags <- metaFrags[ which( metaFrags > 0 ) ]
  
  # check for end of chr
  if(vpIdx > length(frags)-vpSize) {metaFrags <- (length(frags)-vpSize+1):length(frags)}
  metaFrags <- metaFrags[ which( metaFrags <= length(frags) ) ]
  
  if ( length( metaFrags ) != vpSize ){
    # we check if we are at the start or the end of a chromosome
    if ( max(metaFrags) < vpSize ) {
      metaFrags <- seq( 1, vpSize, 1 )
    }
  }
  if ( max(metaFrags) > fragsChr$fragID[length(fragsChr)] ){
    metaFrags <- seq( fragsChr$fragID[length(fragsChr)]-vpSize+1, fragsChr$fragID[length(fragsChr)], 1 )
  } 
  
  # take the viewpoint extremes for computing vp coverage
  vpGR <- reduce(fragsChr[ metaFrags ])
  minVPpos <- min( start(vpGR) ) - 100000
  maxVPpos <- max( end( vpGR )) + 100000
  
  if ( minVPpos < 1 ) {
    minVPpos <- 1
  }
  
  if ( maxVPpos > max( end( fragsChr ) ) ){
    maxVPpos <- max( end( fragsChr ) )
  }

    # selection of the reads that share the viewpoint fragments
  idxPE1.c1 <- queryHits( findOverlaps(reads$PE1,fragsChr[metaFrags]) )
  idxPE2.c1 <- queryHits( findOverlaps(reads$PE2,fragsChr[metaFrags]) )
  
  # the position of the fragment is set to be in the center.
  vpF$pos <- start( resize( vpF, width=1, fix="center") )
  
  # now we sum up all the reads found in each position in range that share one of the ends with the viewpoint fragment
  if ( ( length( vpF ) > 0 ) & ( length( idxPE1.c1 ) > 0 ) & ( length( idxPE2.c1 ) > 0 ) ) {
    vpF$reads <- countOverlaps( vpF,reads$PE2[ idxPE1.c1 ] ) + countOverlaps( vpF,reads$PE1[ idxPE2.c1 ] )
    if ( length( which( vpF$reads > 0 ) ) > 2 ) {
      vpF <- normV4C( vpF )
      vpF$normV4C <- zoo::rollmean( x=vpF$normReads, k=k, fill="extend" )
    } else {
      vpF$normReads <- 0
      vpF$normV4C <- 0
    }
  } else {
    vpF$reads <- 0
    vpF$normReads <- 0
    vpF$normV4C <- 0
  }
  
  vpCovGR <- vpF[ which( start( vpF ) >= minVPpos & end( vpF ) <= maxVPpos ) ]
  #######################why take 100 here? Should be read length?##############################################
  vp_cov <- 100*( length( which( vpCovGR$reads > 0 ) ) / length( vpCovGR ) )
  
  nReads <- length( which( vpF$reads > 0 ) )
  
  if ( nReads != 0 ){
    plotCov <- 100*( nReads / length( vpF ) )
  } else {
    plotCov <- 0
  }
  
  return( list( GR=vpF, vpGR=vpGR, vpCov=vp_cov, plotCov=plotCov ) )
}

normV4C <- function( readsGR, nReads=10e3, nTop=1 ) {
  readsGR$normReads <- 0
  sumTop <- sum( -sort( -readsGR$reads )[ 1:nTop] )
  wNorm <- nReads/( sum( readsGR$reads )-sumTop )
  readsGR$normReads <- wNorm*readsGR$reads
  return(readsGR)
}

doPeakCCall <- function(vpID,peakHiCObj,configOpt,wSize=NULL,alphaFDR=NULL,qWr=NULL,minDist=NULL) {
  
  if(is.null(wSize)){
    
    wSize <- configOpt$peakCalls$wSize
    
  }
  
  if(is.null(alphaFDR)){
    
    alphaFDR <- configOpt$peakCalls$alphaFDR
    
  }
  
  if(is.null(qWr)){
    
    qWr <- configOpt$peakCalls$qWr
    
  }
  
  if(is.null(minDist)){
    
    minDist <- configOpt$peakCalls$minDist
    
  }
  
  rdsFldr <- paste0(configOpt$projectFolder,"rds/")
  hicCond <- configOpt$hicCond
  
  vpsGR <- peakHiCObj[["vpsGR"]]
  partID <- vpsGR$partID[match(vpID,vpsGR$vpID)]
  vpPos <- start(vpsGR[match(vpID,vpsGR$vpID)])
  
  fRDS <- paste0(rdsFldr,"vpReads_",partID,".rds")
  
  if(file.exists(fRDS)) {
    
    vpReads <- readRDS(fRDS)
    nextVPReads <- vpReads[[vpID]][[hicCond]]
    
    num.exp <- length(nextVPReads)
    
    db <- list()
    
    for(i in 1:num.exp) {
      
      if(!is.null(nextVPReads[[i]]$GR)){
        db[[i]] <- data.frame(pos=nextVPReads[[i]]$GR$pos,reads=nextVPReads[[i]]$GR$reads)
      }
      
    }
    if(length(db)==num.exp) {
      resPeakC <- combined.analysis(data=db,num.exp=num.exp,vp.pos=vpPos,wSize=wSize,alphaFDR=alphaFDR,qWr=qWr,minDist=minDist)
      return(resPeakC)
    }
    
  } else {
    
  }
  
}


combine.reps <- function (data) {
  
  num.exp <- length(data)
  
  data.m <- data[[1]]
  
  for (i in 2:num.exp) {
    data.m <- base::merge(data.m, data[[i]], by = 1)
  }
  
  sumReads <- apply(data.m[,2:(num.exp+1)],1,sum)
  out <- data.frame(pos=data.m[,1],reads=sumReads)
  
  return(out)
}

peakAnalysis <- function (data, num.exp = 3, vp.pos, wSize = 21, alphaFDR = 0.1,qWr = 1, minDist = 25e3) {
  #browser()
  if (num.exp == 0) {
    num.exp = data$num.exp
  }
  if (length(vp.pos) == 1) {
    vp.pos <- c(vp.pos, vp.pos)
  }
  vp.pos <- sort(vp.pos)
  db <- tryCatch({ combine.experiments(data, num.exp, vp.pos) }, error=function(e) { return (NULL) } ) #peakC::combine.experiments   copy-pasted now
  if ( !is.null( db ) ){
    dbR <- db
    dbR[, 2:(num.exp + 1)] <- apply(db[, 2:(num.exp + 1)], 2, 
                                    caTools::runmean, k = wSize, endrule = "mean")
    dbR[, 2:(num.exp + 1) + num.exp] <- apply(db[, 2:(num.exp + 
                                                        1) + num.exp], 2, caTools::runmean, k = 5, endrule = "mean") #why k=5 here?

    pseudoCount <- apply(db[, 2:(num.exp + 1)], 2, non.zero.quantile, 
                         probs = 0.05)
    # print(pseudoCount)
    if(!is.na(sum(pseudoCount) ) ){
      pseudoCount <- sum(pseudoCount)/num.exp
      ratio <- cbind(db[, 1], (dbR[, 2:(num.exp + 1)] + pseudoCount)/(dbR[, 
                                                                          (2:(num.exp + 1)) + num.exp] + pseudoCount))
      delta <- cbind(db[, 1], dbR[, 2:(num.exp + 1)] - dbR[, (2:(num.exp + 
                                                                   1)) + num.exp])
      p.val <- rank.product.p(data = dbR, num.exp = num.exp, method = "diff")
      sfr <- significant.fragments(p.value = p.val, pos = db[, 
                                                             1], window = wSize, FDR = alphaFDR)
      sel.frag <- db[which((db[, 1] < vp.pos[1] & vp.pos[1] - db[, 
                                                                 1] > minDist) | (db[, 1] > vp.pos[2] & db[, 1] - vp.pos[2] > 
                                                                                    minDist)), 1]
      idx <- delta[, 1] %in% sel.frag
      tfr <- thresholdFrags(resids = apply(ratio[idx, 2:(num.exp + 
                                                           1)], 1, mean), frags = ratio[idx, 1], wSize = wSize, 
                            qW = qWr)
      sfr <- intersect(sfr, tfr)
      list(dbR = dbR, peak = sfr, num.exp = num.exp, p.value = p.val, 
           ratio = apply(ratio[, 2:(num.exp + 1)], 1, mean), delta = apply(delta[, 2:(num.exp + 1)], 1, mean), sel = sel.frag)
    } else {
      return(NULL)
    }
  }
}

significant.fragments <- function( p.value, pos, window = 21, FDR = 0.01 ){
  #correct the nominal p-value for multiple hypothesis testing
  p.combined <- p.adjust(p.value, method="fdr")
  #determine the significant windows and select the fragments therein
  sig.i <- which(p.combined < FDR)
  if(length(sig.i)>0) {
    sig.i.start <- sig.i-floor(window/2); sig.i.end <- sig.i+floor(window/2)
    sig.i <- unique(multi.seq(sig.i.start,sig.i.end))
    sig.i <- sig.i[sig.i >= 1 & sig.i <= length(pos)]
    sigFrags <- pos[sig.i]
    return(sigFrags)
  } else {
    return(NULL)
  }
}

righttailgamma <- function(r,k,n) {
  
  return(1 - pgamma(-log(r/(n+1)^k),k,scale=1))
  
}

rank.product.p <- function( data, num.exp,method="diff"){
  if(method=="diff") {
    stats <- data[,2:(num.exp+1)]-data[,(2:(num.exp+1))+num.exp]
  } else {
    stats <- data[,2:(num.exp+1)]/data[,(2:(num.exp+1))+num.exp]
  }
  rp <- nrow(data)-apply(stats,2,rank)+1
  rp <- apply(rp,1,prod)
  p <- righttailgamma(rp,num.exp,length(rp))
}

thresholdFrags <- function(resids,frags,wSize=21,qW=5) {
  
  qMax <- getThreshold(resids=resids,qW=qW)
  
  sel.i <- which(resids > qMax)
  if(length(sel.i)>0) {
    sel.i.start <- sel.i-floor(wSize/2); sel.i.end <- sel.i+floor(wSize/2)
    sel.i <- unique(multi.seq(sel.i.start,sel.i.end))
    sel.i <- sel.i[sel.i >= 1 & sel.i <= length(frags)]
    selFrags <- frags[sel.i]
    return(selFrags)
  }else{
    return(NULL)
  }
  
  
}

non.zero.quantile <- function( x, probs ){
  if(sum(is.na(x))==0) {
    
    return(quantile(x[x > 0], probs))
    
  } else {
    
    return(NA)
    
  }
}

rem <- function(a, n ){
  half.window <- floor(n/2)
  head(tail(a, -half.window),-half.window)
}

#quick way of generating a vector with the required indexes
multi.seq <- function( start, end ){
  x <- rep(start, end-start+1)->x
  df <- diff(x)
  df <- df + 1
  low <- which(df > 1)
  df[low] <- -diff(c(0,low))+1
  add <- c(0,cumsum(df))
  x + add
}

getThreshold <- function(resids,qW=5) {
  q75 <- quantile(resids,probs=0.75) #75% quantile of the residuals
  qd50 <- diff(quantile(resids,probs=c(0.25,0.75))) #the range between the 25% and 75% quantiles
  threshold <- q75 + qW*qd50
  return(threshold)
}

collapseFrags <- function(frags, peakFrags) {
  
  N <- length(frags)
  
  if(length(peakFrags)>0){
    
    fragRanges <- IRanges(frags,frags)
    end(fragRanges)[1:(N-1)] <- start(fragRanges)[2:N]-1
    
    binaryFrags <- ifelse(frags%in%peakFrags,1,0)
    
    return(reduce(split(fragRanges,binaryFrags)[["1"]],min.gapwidth=1))
    
  } else {
    
    return(NULL)
  }
  
}

getPartitionPeaks <- function(partID, peakHiCObj, configOpt=NULL, hicCond=NULL, wSize=NULL, qWr=NULL, alphaFDR=NULL,  minDist=NULL, nReps=NULL, writePeaksFile=FALSE) {
  #browser()
  if(is.null(configOpt)){
    configOpt <- peakHiCObj$configOpt
  }
  
  if(is.null(hicCond)){
    hicCond <- configOpt$hicCond
  }
  
  if(is.null(nReps)) {
    nReps <- sum(as.vector(peakHiCObj[["hic"]][["design"]][["HiCMap"]])==hicCond)
  }
  
  if(is.null(wSize)){
    wSize <- configOpt$peakCalls$wSize
  }
  
  if(is.null(qWr)){
    qWr <- configOpt$peakCalls$qWr
  }
  
  if(is.null(alphaFDR)){
    alphaFDR <- configOpt$peakCalls$alphaFDR
  }
  
  if(is.null(minDist)){
    minDist <- configOpt$peakCalls$minDist
  }
  
  rdsFldr <- paste0(configOpt$projectFolder,"rds/")
  profilesFldr <- paste0(rdsFldr,"profiles/")
  loopsFldr <- paste0(rdsFldr,"loops/")
  
  vpsGR <- peakHiCObj[["vpsGR"]]
  partVPs <- vpsGR[vpsGR$partID==partID]
  
  if(length(partVPs)>0) {
    
    peakHeader <- paste("chr","vp_X1","vp_X2","anchor_X1","anchor_X2","anchorSize","maxV4Cscore","maxV4CscorePos","delta","ratio","minPval","id",sep="\t")
    
    if(writePeaksFile){
      
      loopFile <- paste0(loopsFldr,peakHiCObj$name,"_GW_nReps_",nReps,"_peakHiC_wSize_",wSize,"_qWr_",qWr,"_alphaFDR_",alphaFDR,"_loops.txt")
      
    } else {
      
      loopFile <- tempfile()
      
    }
    
    if(!file.exists(loopFile)) {
      
      cat(c(peakHeader,"\n"),file=loopFile,sep="")
      
    }
    
    profilesFile <- paste0(profilesFldr,partID,".rds")
    
    vpChr <- as.vector(seqnames(vpsGR[match(partID,vpsGR$partID)]))[1]
    
    if(file.exists(profilesFile)) {
      
      V4Cs <- readRDS(profilesFile)
      
    } else {
      
      V4Cs <- list()
      
    }
    
    V4Cs[[hicCond]] <- list()
    
    
    if(length(partVPs)>0) {
      
      fRDS <- paste0(rdsFldr,"vpReads_",partID,".rds")
      ids <- partVPs$vpID
      
      if(file.exists(fRDS)) {
        
        vpReads <- readRDS(fRDS)
        
        #only executes this part if there are no ids found in the partitioned V4C read files, e.g. in vpReads_part1.rds
        if(sum(!ids%in%names(vpReads))>0){
          
          fragsFile <- configOpt$fragsFile
          frags <- readRDS(fragsFile)[[vpChr]]
          vpReads <- getPeakHiCData(partID=partID,frags=frags,peakHiCObj=peakHiCObj)
          
        }
        
      } else { #only executes this if no partitioned V4C files were meade (i.e. callPartitionPeaksbyChr.R, e.g. vpReads_part1.rds)
        
        fragsFile <- configOpt$fragsFile
        frags <- readRDS(fragsFile)[[vpChr]]
        vpReads <- getPeakHiCData(partID=partID,frags=frags,peakHiCObj=peakHiCObj)
        
      }
      
      for(id in ids) {
        #consider subsetting vpReads before passing to function
        peakRes <- getPeakCPeaksWithReps(vpID=id,vpReads=vpReads, partVPs=partVPs, hicCond=hicCond, wSize=wSize, qWr=qWr, alphaFDR=alphaFDR,  minDist=minDist)
        
        if(nrow(peakRes$df)>0) {
          
          write.table(peakRes$df,file=loopFile,append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
          
        }
        
        V4Cs[[hicCond]][[id]] <- peakRes$gR
        
      }
      
      saveRDS(V4Cs,file=profilesFile)
      
    }
    
    if(!writePeaksFile) {
      
      peakRes <- read.table(loopFile,sep="\t",stringsAsFactors=FALSE,header=TRUE)
      return(peakRes)
      
    }
  }
}

# getChrPeaks <- function(chr, peakHiCObj, configOpt, hicCond=NULL, wSize=NULL, qWr=NULL, alphaFDR=NULL,  minDist=NULL, nReps=NULL, nThreads=8) {
#   #do we need doParallel can we even do this on top of exisiting clusters?
#   suppressPackageStartupMessages(require(doParallel))
#   registerDoParallel(cores=nThreads)
#   
#   vpsGR <- peakHiCObj[["vpsGR"]]
#   ids <- unique(subChr(vpsGR,chr)$partID)
#   
#   foreach(i = 1:length(ids)) %dopar% {
#  # for(i in 1:length(ids)){
#     getPartitionPeaks(partID=ids[i], peakHiCObj=peakHiCObj, hicCond=hicCond, wSize=wSize, qWr=qWr, alphaFDR=alphaFDR,  minDist=minDist, nReps=nReps)
#     
#   }
#   
#   stopCluster()
#   
# }

getPeakCPeaksWithReps <- function( vpID, vpReads, partVPs, hicCond="HAP1_WAPL", wSize=31, qWr=1, alphaFDR=0.1,  minDist=30e3 ) {

  vpPos <- start(partVPs[match(vpID,partVPs$vpID)])
  vpChr <- as.vector(seqnames(partVPs[match(vpID,partVPs$vpID)]))
  
  nextVPReads <- vpReads[[vpID]][[1]]
  
  num.exp <- length(nextVPReads)
  
  peakCReads <- list()
  
  for(i in 1:num.exp) {
    
    if(!is.null(nextVPReads[[i]]$GR)){
      dat <- data.frame(pos=nextVPReads[[i]]$GR$pos,reads=nextVPReads[[i]]$GR$reads)
      colnames(dat) <- c("pos",paste0("reads.R",i))
      peakCReads[[i]] <- dat
      
    } else {
      
      peakGR <- NULL
      maxV4C <- NULL
      minPval <- NULL
      maxRatio <- NULL
      maxDelta <- NULL
      finalGR <- NULL
      maxV4CscorePos <- NULL
      finalDF <- matrix(0,nrow=0,ncol=1)
      
      return( list( gR=finalGR, df=finalDF ) )
    }
  }
  #browser()
  peakCReads$num.exp <- num.exp
  nReps <- num.exp
  
  peakCRes <- peakAnalysis( data=peakCReads, num.exp=nReps, vp.pos=vpPos, wSize=wSize, minDist=minDist, qWr=qWr, alphaFDR=alphaFDR )
  
  if( length( peakCRes$peak ) > 0 ) {
    
    peakGR <- sort(collapseFrags(frags=peakCRes$dbR[,1],peakFrags=peakCRes$peak))
    y.ave <- apply(peakCRes$dbR[, 2:(nReps + 1)], 1, median)
    gR <- IRanges(peakCRes$dbR$pos,peakCRes$dbR$pos)
    finalGR <- GRanges( unique( vpChr ), gR, normV4C=y.ave )
    ovl <- findOverlaps(gR,peakGR)
    maxV4C <- as.vector(tapply(y.ave[ovl@from],ovl@to,max))
    minPval <- as.vector(tapply(peakCRes$p.value[ovl@from],ovl@to,min))
    maxRatio <- as.vector(tapply(peakCRes$ratio[ovl@from],ovl@to,max))
    maxDelta <- as.vector(tapply(peakCRes$delta[ovl@from],ovl@to,max))
    ttt <- data.frame( gR, y.ave )[queryHits(ovl),]
    ttt$loop <- subjectHits( ovl )
    ttt <- split( ttt, ttt$loop )
    maxV4CscorePos <- sapply( 1:length( ttt ), function(x) ttt[[x]]$start[ which.max( ttt[[x]]$y.ave ) ] )
    finalDF <- data.frame( 
      chr=vpChr
      , vp_X1=vpPos
      , vp_X2=vpPos
      , anchor_X1=start( peakGR )
      , anchor_X2=end( peakGR )
      , anchorSize=width( peakGR )
      , maxV4Cscore=maxV4C
      , maxV4CscorePos=maxV4CscorePos
      , delta=maxDelta
      , ratio=maxRatio
      , minPval=minPval
      , id=vpID
      , stringsAsFactors=FALSE
    )
  } else {
    
    peakGR <- NULL
    maxV4C <- NULL
    minPval <- NULL
    maxRatio <- NULL
    maxDelta <- NULL
    finalGR <- NULL
    maxV4CscorePos <- NULL
    finalDF <- matrix(0,nrow=0,ncol=1)
    
  }  
  
  return( list( gR=finalGR, df=finalDF ) )
  
}

doPeakCPlot <- function(vpID,peakHiCObj,configOpt,hicCond=NULL,wSize=NULL,qWr=NULL,alphaFDR=NULL,minDist=NULL) {
  
  if(is.null(hicCond)){
    hicCond <- configOpt$hicCond
  }
  
  if(is.null(wSize)){
    wSize <- configOpt$peakCalls$wSize
  }
  
  if(is.null(qWr)){
    qWr <- configOpt$peakCalls$qWr
  }
  
  if(is.null(alphaFDR)){
    alphaFDR <- configOpt$peakCalls$alphaFDR
  }
  
  if(is.null(minDist)){
    minDist <- configOpt$peakCalls$minDist
  }
  
  rdsFldr <- paste0(configOpt$projectFolder,"rds/")
  
  vpsGR <- peakHiCObj[["vpsGR"]]
  partID <- vpsGR$partID[match(vpID,vpsGR$vpID)]
  vpPos <- start(vpsGR[match(vpID,vpsGR$vpID)])
  
  fRDS <- paste0(rdsFldr,"vpReads_",partID,".rds")
  
  if(file.exists(fRDS)) {
    
    vpReads <- readRDS(fRDS)
    nextVPReads <- vpReads[[vpID]][[hicCond]]
    
  } else {
    
  }
  
  num.exp <- length(nextVPReads)
  
  db <- list()
  db$num.exp <- num.exp
  
  for(i in 1:num.exp) {
    
    db[[i]] <- data.frame(pos=nextVPReads[[i]]$GR$pos,reads=nextVPReads[[i]]$GR$reads)
    
  }
  
  resPeakC <- combined.analysis(data=db,num.exp=num.exp,vp.pos=vpPos,wSize=wSize,alphaFDR=alphaFDR,qWr=qWr,minDist=minDist)
  
  return(resPeakC)
  
}

getPeakHiCData <- function(partID,frags,peakHiCObj,configOpt=NULL,hicCond=NULL,wSize=NULL,vpSize=NULL,v4cSize=NULL) {
  #browser()
  if(is.null(configOpt)){
    
    configOpt <- peakHiCObj$configOpt
    
  }
  
  if(is.null(hicCond)){
    hicCond <- configOpt$hicCond
  }
  
  if(is.null(wSize)){
    wSize <- configOpt$V4C$wSize
  }	
  
  if(is.null(vpSize)){
    vpSize <- configOpt$V4C$vpSize
  }
  
  if(is.null(v4cSize)){
    v4cSize <- configOpt$V4C$v4cSize
  }
  
  designMat <- peakHiCObj[["hic"]][["design"]]
  hicTracksByCondition <- split(as.vector(designMat$trackID),designMat$HiCMap)
  
  vps <- getVPs(partID,peakHiCObj)
  vpReads <- list()
  
  if(length(vps$vpID)>0) {
    
    reads <- list()
    tracks <- hicTracksByCondition[[hicCond]]
    if(is.null(tracks)){message("peakHiC ERROR: hicCond from config.yml does not coincide with any condition from design file")
      q('no')}
    
    #browser()
    for(trackID in tracks) {
      
      reads[[trackID]] <- getPartitionReads(partID=partID,trackID=trackID,peakHiCObj=peakHiCObj)
      
    }
    
    for(vpID in vps$vpID) {
     
      vpReads[[vpID]] <- list()
      vpReads[[vpID]][[hicCond]] <- list()
      fragID <- vps$fragID[match(vpID,vps$vpID)]
      vpPos <- start(vps[match(vpID,vps$vpID)])
    
    for(trackID in tracks) {
     
      vpReads[[vpID]][[hicCond]][[trackID]] <- extractReads(fragID=fragID,vpS=vpPos,vpZoom=c(v4cSize,v4cSize),reads=reads[[trackID]],frags=frags,k=wSize,vpSize=vpSize)
     
    }
    }

  }
  
  vpReads$configOpt <- configOpt
  return(vpReads)
  
}

vpID.V4C <- function(vpID, vps, hicCond, tracks){
  #per viewpoint contruct a V4C
  #browser()
  sub.vpReads <- list()
  sub.vpReads[[hicCond]] <- list()
  fragID <- vps$fragID[match(vpID,vps$vpID)]
  vpPos <- start(vps[match(vpID,vps$vpID)])
  
  #ptm <- proc.time()
  #per track (per VP) gather the reads
  sub.vpReads[[hicCond]] <- lapply(tracks, FUN=function(trackID) extractReads(
    fragID=fragID,vpS=vpPos,vpZoom=c(v4cSize,v4cSize),reads=reads[[trackID]],frags=frags,k=wSize,vpSize=vpSize) )
  #print(proc.time()-ptm)
  
  names(sub.vpReads[[hicCond]]) <- tracks
  
  return(sub.vpReads)
}  


subChr <- function(gR,chr) {
  
  idx <- as.vector(seqnames(gR))==chr
  return(gR[idx])
  
}

readHiCDesign <- function(designFile) {
  
  designDF <- read.table(designFile,stringsAsFactors=FALSE,header=TRUE)
  cols <- colnames(designDF)
  
  if(length(cols)==3) {
    
    if(sum(colnames(designDF)==c("HiCMap","repID","trackID"))){
      
      return(designDF)  
    
    } else {
      
      message("HiC design matrix should contain 3 columns named HiCMap, repID and trackID")
      return(NULL)
      
    }
    
  } else {
    
    message("HiC design matrix should contain 3 columns named HiCMap, repID and trackID")
    return(NULL)
    
  }
  
}

cpet <- function(i,Sq1,Sq2){return(length(intersect(Sq1[[i]],Sq2[[i]])))}

tagRegions <- function(PE1,PE2,qR1,qR2) {
  
  olapPE1q1 <- findOverlaps(qR1,PE1)
  olapPE2q2 <- findOverlaps(qR2,PE2)
  
  olapPE2q1 <- findOverlaps(qR1,PE2)
  olapPE1q2 <- findOverlaps(qR2,PE1)
  
  Sq1 <- split(olapPE1q1@to,factor(olapPE1q1@from,1:length(qR1)))
  Sq2 <- split(olapPE2q2@to,factor(olapPE2q2@from,1:length(qR2)))
  
  tags <- sapply(1:length(qR1),cpet,Sq1=Sq1,Sq2=Sq2)
  
  Sq1 <- split(olapPE1q2@to,factor(olapPE1q2@from,1:length(qR1)))
  Sq2 <- split(olapPE2q1@to,factor(olapPE2q1@from,1:length(qR2)))
  
  tags <- tags+sapply(1:length(qR1),cpet,Sq1=Sq1,Sq2=Sq2)
  
  return(tags)
  
}

getLoopCovbyChr <- function(chr,loops,peakHiCObj,configOpt=NULL,hicCond=NULL,anchorSize=10e3,loopCovFile=NULL,writePeaksFile=FALSE) {
  
  vpsGR <- peakHiCObj$vpsGR
  vpIDs <- loops$id[loops$chr==chr]
  
  if(length(vpIDs)>0) {
    
    partIDs <- unique(vpsGR$partID[match(vpIDs,vpsGR$vpID)])
    
    if(length(partIDs)>0) {
      
      for(id in partIDs) {
        
        getLoopCovbyPartition(partID=id,loops=loops,peakHiCObj=peakHiCObj,configOpt=configOpt,hicCond=hicCond,anchorSize=anchorSize,loopCovFile=loopCovFile,writePeaksFile=writePeaksFile)
        
      }
    }
    
  }
  
}

getLoopCovbyPartition <- function(partID,loops,peakHiCObj,configOpt=NULL,hicCond=NULL,anchorSize=10e3,loopCovFile=NULL,writePeaksFile=FALSE) {
  
  if(is.null(configOpt)){
    
    configOpt <- peakHiCObj$configOpt
    
  }
  
  if(is.null(hicCond)){
    
    hicCond <- configOpt$hicCond
  }
  
  wSize <- configOpt$peakCalls$wSize
  qWr <- configOpt$peakCalls$qWr
  alphaFDR <- configOpt$peakCalls$alphaFDR
  
  cNames <- colnames(loops)
  if(sum(!c("loopID","id")%in%cNames)){stop(paste0("loop data.frame must have loopID (identifier of loop) and id (identifier of viewpoint) columns"))}
  
  designMat <- peakHiCObj[["hic"]][["design"]]
  hicTracksByCondition <- split(as.vector(designMat$trackID),designMat$HiCMap)
  tracks <- hicTracksByCondition[[hicCond]]
  
  nReps <- length(tracks)
  
  rdsFldr <- paste0(configOpt$projectFolder,"rds/")
  loopsFldr <- paste0(rdsFldr,"loops/")
  
  if(is.null(loopCovFile)) {
    
    loopCovFile <- paste0(loopsFldr,peakHiCObj$name,"_GW_nReps_",nReps,"_peakHiC_wSize_",wSize,"_qWr_",qWr,"_alphaFDR_",alphaFDR,"_loopCov.txt")
    
  }
  
  vps <- getVPs(partID,peakHiCObj)$vpID
  
  if(length(vps)>0) {
    
    partLoops <- loops[loops$id%in%vps,]
    
    if(nrow(partLoops)>0) {
      
      lVP <- floor((partLoops$vp_X1+partLoops$vp_X2)/2)
      lAnchor <- partLoops$maxV4CscorePos
      xPos <- pmin(lVP,lAnchor)
      yPos <- pmax(lVP,lAnchor)
      lChr <- partLoops$chr
      lx <- resize(GRanges(lChr,IRanges(xPos,xPos)),width=anchorSize,fix="center")
      ly <- resize(GRanges(lChr,IRanges(yPos,yPos)),width=anchorSize,fix="center")
      
      out <- as.data.frame(matrix(0,nrow=length(lx),ncol=(length(tracks)+1)))
      colnames(out) <- c("loopID",tracks)
      out$loopID <- partLoops$loopID
      
      for(trackID in tracks) {
        
        PRreads <- getPartitionReads(partID=partID,trackID=trackID,peakHiCObj=peakHiCObj)
        loopTags <- tagRegions(PE1=PRreads$PE1,PE2=PRreads$PE2,qR1=lx,qR2=ly)
        out[[trackID]] <- loopTags
        
      }
      
      loopsOut <- cbind(loops[match(out$loopID,loops$loopID),],out[,-1])
      
      if(writePeaksFile) {
        
        fileHeader <- paste0(colnames(loopsOut),collapse="\t")
        
        if(!file.exists(loopCovFile)){cat(c(fileHeader,"\n"),file=loopCovFile,sep="")}
        
        write.table(loopsOut,file=loopCovFile,append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
        
      } else{
        
        return(loopsOut)
        
      }
    }
  }
  
}

normalizeLoopCov <- function(loops,hicCond="peakHiC",nBins=150) {
  
  covVar <- paste0(hicCond,".covQ")
  covNormVar <- paste0(hicCond,".covQ.norm")
  vp0 <- floor((loops$vp_X1+loops$vp_X2+1)/2)
  
  loops$dist <- abs(vp0-loops$maxV4CscorePos)
  loops <- loops[order(loops$dist,decreasing=TRUE),]
  
  distBins <- rep(1:nBins,each=ceiling(nrow(loops)/nBins))[1:nrow(loops)]
  covMean <- tapply(loops[[covVar]],distBins,mean)
  loops[[covNormVar]] <- loops[[covVar]]-covMean[distBins]
  
  return(loops)
  
}

binGenome <- function(BSgenome,binSize=10e3){
  
  do.call( require, args=list( BSgenome ) )
  assign( 'genome', base::get(as.vector(BSgenome)))
  
  x<-seqinfo(genome)
  bins <- tileGenome(x, tilewidth=binSize, cut.last.tile.in.chrom=TRUE)
  
  return(bins)
  
}

getBinPairs <- function(loops,bins,anchorSize=10e3) {
  
  lx <- resize(GRanges(seqnames=loops$chr,IRanges(loops$vp_X1,loops$vp_X2)),width=1,fix="center")
  ly <- resize(GRanges(seqnames=loops$chr,IRanges(loops$maxV4CscorePos,loops$maxV4CscorePos)),width=1,fix="center")
  
  ovl.lx <- findOverlaps(lx,bins)
  ovl.ly <- findOverlaps(ly,bins)
  
  loopIDs <- unique(intersect(ovl.lx@from,ovl.ly@from))
  loops <- loops[loopIDs,]
  
  ovl.lx <- ovl.lx[match(loopIDs,ovl.lx@from)]
  ovl.ly <- ovl.ly[match(loopIDs,ovl.ly@from)]
  
  NR.binID <- paste0(pmin(ovl.lx@to,ovl.ly@to),"_",pmax(ovl.lx@to,ovl.ly@to))
  loops$NR.binID <- NR.binID
  
  return(loops)
  
}


getTADCovbyPartition <- function(partID,loops,peakHiCObj,hicCond="Rao_4DN_GM12878",anchorSize=10e3) {
  
  if(is.null(loops$loopID)|is.null(loops$id)){stop(paste0("loop data.frame must have loopID (identifier of loop) and id (identifier of viewpoint) columns"))}
  
  designMat <- peakHiCObj[["hic"]][["design"]]
  hicTracksByCondition <- split(as.vector(designMat$trackID),designMat$HiCMap)
  
  vps <- getVPs(partID,peakHiCObj)$vpID
  partLoops <- loops[loops$id%in%vps,]
  
  if(nrow(partLoops)>0) {
    
    tracks <- hicTracksByCondition[[hicCond]]
    lVP <- floor((partLoops$vp_X1+partLoops$vp_X2)/2)
    lAnchor <- partLoops$maxV4CscorePos
    xPos <- pmin(lVP,lAnchor)
    yPos <- pmax(lVP,lAnchor)
    lChr <- partLoops$chr
    lx <- resize(GRanges(lChr,IRanges(xPos,xPos)),width=anchorSize,fix="center")
    ly <- resize(GRanges(lChr,IRanges(yPos,yPos)),width=anchorSize,fix="center")
    
    out <- as.data.frame(matrix(0,nrow=length(lx),ncol=(length(tracks)+1)))
    colnames(out) <- c("loopID",tracks)
    out$loopID <- partLoops$loopID
    
    for(trackID in tracks) {
      
      PRreads <- getPartitionReads(partID=partID,trackID=trackID,peakHiCObj=peakHiCObj)
      loopTags <- tagRegions(PE1=PRreads$PE1,PE2=PRreads$PE2,qR1=lx,qR2=ly)
      out[[trackID]] <- loopTags
      
    }
  }
  
  
  return(out)
  
}

getPartitionReads <- function(partID,trackID,peakHiCObj,configOpt=NULL,nThread=2) {
  require(data.table)
  require(GenomicRanges)
  
  if(is.null(configOpt)) {
        configOpt <- peakHiCObj$configOpt
    }
  
  readsFldr <- configOpt$hicReadsFldr
  pairixBinary <- configOpt$pairixBinary
  
  out <- NULL
  
  partIdx <- match(partID,peakHiCObj$partition$partGR$partID)
  
  if(!is.na(partIdx)){
    
    qPartition <- peakHiCObj$partition$partGR[partIdx]
    chr <- as.vector(seqnames(qPartition[1]))
    qRegion <- paste0("\'",chr,":",start(qPartition),"-",end(qPartition),"\'")
    tmpFldr <- tempdir()
    tmpFile <- paste0(tmpFldr,"/pairix_query.txt")
    pairixFile <- paste0(readsFldr,'/',trackID,".pairs.gz")
    #browser()
    if (file.exists(pairixFile)&file.exists(pairixBinary)) {
      
      dir.create(tempdir(), showWarnings = FALSE)
      
      cmd <- paste0(pairixBinary," ",pairixFile," ",qRegion," > ",tmpFile)
      system(cmd)
      
      dat <- fread(file=tmpFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,nThread=nThread)
      
      cmd <- paste0("rm -r ",tempdir())
      system(cmd)

      #Now we save reads in a dedicated location and remove right after returning GRanges object
      #we could use R package of pairix thereby circumventing file creation altoghether:
      #query per VP the pairix file to obtain reads as opposed to using intermediate parition GR file
      PE1 <- GRanges(seqnames=dat$V2,IRanges(dat$V3,dat$V3))
      PE2 <- GRanges(seqnames=dat$V4,IRanges(dat$V5,dat$V5))
      out <- list(PE1=PE1,PE2=PE2)
      
    } else{message('\npeakHiC WARNING: pairix location ',pairixBinary,' or reads file ',pairixFile,' not found.\n')}
    
  } else{message('\npeakHiC WARNING: partition ',partID, ' not mentioned in peakHiCObj.\n')}
  
  return(out)
  
}

getVPs <- function(partID,peakHiCObj) {
  
  if(is.null(peakHiCObj[["vpsGR"]])) {
    
    out <- NULL
    
  } else {
    
    vpsGR <- peakHiCObj[["vpsGR"]]
    out <- vpsGR[vpsGR$partID==partID]
    
  }
  
  return(out)
}

getV4CData <- function(vpID,peakHiCObj,configOpt=NULL,wSize=21,alphaFDR=0.1,qWr=2.0,minDist=20e3,storeVPReads=FALSE) {
  
  if(is.null(configOpt)) {
    
    configOpt <- peakHiCObj$configOpt
    
  }
  
  rdsFldr <- paste0(configOpt$projectFolder,"rds/")
  fragsFile <- configOpt$fragsFile
  hicCond <- configOpt$hicCond
  
  vpsGR <- peakHiCObj[["vpsGR"]]
  partID <- vpsGR$partID[match(vpID,vpsGR$vpID)]
  vpPos <- start(vpsGR[match(vpID,vpsGR$vpID)])
  vpChr <- as.vector(seqnames(vpsGR[match(vpID,vpsGR$vpID)]))[1]
  
  fRDS <- paste0(rdsFldr,"vpReads_",partID,".rds")
  
  if(file.exists(fRDS)) {
    
    vpReads <- readRDS(fRDS)
    
  } else {
    
    frags <- tryCatch(readRDS(fragsFile)[[vpChr]], error=function(cond){message("\npeakHiC ERROR: Wrong RE fragment file format: fragments should be a I/GRanges object in a list per chromsome.\n")})
    vpReads <- getPeakHiCData(partID=partID,frags=frags,peakHiCObj=peakHiCObj)
    
    if(storeVPReads) {
      
      saveRDS(vpReads,file=fRDS)
      
    }
  }
  
  nextVPReads <- vpReads[[vpID]][[hicCond]]
  
  num.exp <- length(nextVPReads)
  
  db <- data.frame(pos=nextVPReads[[1]]$GR$pos,stringsAsFactors=FALSE)
  tracks <- names(nextVPReads)
  
  for(i in 1:length(tracks)) {
    
    db[[tracks[i]]] <- nextVPReads[[i]]$GR$normV4C
    
  }
  
  db[["V4C"]] <- apply(db[,tracks],1,median)
  
  return(db)
} 

getOvLoops <- function(vpIDs,vpsGR,loops) {
  
  require(GenomicRanges)
  
  out <- NULL
  
  lx <- resize(GRanges(seqnames=loops$chr,IRanges(loops$vp_X1,loops$vp_X2)),width=10e3,fix="center")
  ly <- resize(GRanges(seqnames=loops$chr,IRanges(loops$maxV4CscorePos,loops$maxV4CscorePos)),width=10e3,fix="center")
  
  gRs <- vpsGR[na.omit(match(vpIDs,vpsGR$vpID))]
  
  IDX <- which(countOverlaps(lx,gRs)>0|countOverlaps(ly,gRs)>0)
  
  if(length(IDX)>0){
    
    out <- list(lx=lx[IDX],ly=ly[IDX])
    
  } 
  
  return(out)
  
}

createConfig <- function(confFile) {
  
  if( !suppressMessages(require( "config", character.only=TRUE ) ) ) stop( "Package not found: config" )
  configOpt <- suppressWarnings(config::get(file=confFile))
  configOpt$confFile <- attr(configOpt,"file")
  
  return(configOpt)
  
}

getChrs <- function(BSname) {
  
  do.call(require,args=list(BSname))
  assign( 'BScurrent', base::get( BSname ) )
  allSeqs <- seqnames(BScurrent)
  
  Un <- grep("Un",allSeqs)
  Alt <- grep("_alt",allSeqs)
  Rand <- grep("_random",allSeqs)
  Mit <- grep("chrM",allSeqs)
  
  chrs <- allSeqs[setdiff(1:length(allSeqs),c(Un,Alt,Rand,Mit))]
  
  return(chrs)
  
}

getChrSeqlengths <- function(BSname) {
  
  do.call(require,args=list(BSname))
  assign( 'BScurrent', base::get( BSname ) )
  allSeqs <- seqnames(BScurrent)
  
  Un <- grep("Un",allSeqs)
  Alt <- grep("_alt",allSeqs)
  Rand <- grep("_random",allSeqs)
  Mit <- grep("chrM",allSeqs)
  
  chrs <- allSeqs[setdiff(1:length(allSeqs),c(Un,Alt,Rand,Mit))]
  
  return(seqlengths(BScurrent)[chrs])

}

createVPs <- function(vpData, peakHiCObj, fragsFile=NULL) {
  #vpData <- read.table("DATA/example_data/hg38_4DN_Rao_GM12878_peakHiC_example_VPs.txt",header = TRUE, stringsAsFactors = FALSE)
  if(is.null(fragsFile)) {
    
    fragsFile <- peakHiCObj$configOpt$fragsFile
    
  }
  
  requiredNames <- c("chr","vpPos","vpID","name","type")
  
  if( length(intersect(requiredNames,colnames(vpData)))==length(requiredNames) ) {
    
    vpsGR <- GRanges(seqnames=vpData$chr,IRanges(vpData$vpPos,vpData$vpPos))
    vpsGR$vpID <- vpData$vpID
    vpsGR$name <- vpData$name
    vpsGR$type <- vpData$type
    
    partGR <- peakHiCObj$partition$partGR
    
    to.keep <- unique(findOverlaps(ranges(vpsGR), ranges(partGR))@from)
    thrown <- length(vpsGR[-to.keep]$vpID)
    if(thrown>1){message(paste0('Discarded ',thrown,' VPs outside defined paritions.'))}
    
    vpsGR <- vpsGR[to.keep]
    
    distMap <- distanceToNearest(vpsGR,resize(partGR,width=1,fix="center"))
    
    vpsGR$partID <- "no.part"
    vpsGR$partID[distMap@from] <- partGR$partID[distMap@to]
    vpsGR <- sort(vpsGR[vpsGR$partID!="no.part"])

    fragsList <- readRDS(fragsFile)
    chrs <- intersect(names(fragsList),as.vector(seqnames(vpsGR)))
    if(length(chrs)==0){message(paste0("\npeakHiC ERROR: Fragments File is incorrectly formatted (should be a list of GRanges per chromosome)",
    " or it does not contain any of the chromosomes specified in the viewpoints (VP) file.\n"))}
      
    nextChr <- chrs[1]
    nextFrags <- fragsList[[nextChr]]
    nextVPs <- subChr(gR=vpsGR,chr=nextChr)
    nextVPs$fragID <- (-1)
    
    fragOvl <- findOverlaps(nextVPs,nextFrags)
    nextVPs$fragID[fragOvl@from] <- nextFrags$fragID[fragOvl@to]
    out <- nextVPs
    
    if(length(chrs)>1) {
      
      for(nextChr in chrs[-1]) {
        
        nextFrags <- fragsList[[nextChr]]
        nextVPs <- subChr(gR=vpsGR,chr=nextChr)
        nextVPs$fragID <- (-1)
        
        fragOvl <- findOverlaps(nextVPs,nextFrags)
        nextVPs$fragID[fragOvl@from] <- nextFrags$fragID[fragOvl@to]
        out <- c(out,nextVPs)
      }
    }
    
    return(out)
    
  } else {
    
    message("Required colnames in vpData are : chr,vpPos,vpID,name,type ")
    message("returning NULL .....")
    
    return(NULL)
    
  }
  
}

createFrags <- function(configOpt,writeToRDS=FALSE) {
  
  frags <- list()
  
  BSname <- strsplit(configOpt$BSgenome,split=" ")[[1]][[2]]
  
  do.call(require,args=list(BSname))
  assign( 'BScurrent', base::get( BSname ) )
  
  REmotif <- strsplit(configOpt$RE,split=" ")[[1]][[2]]

  chrs <- getChrs(BSname=BSname)
  sLengths <- getChrSeqlengths(BSname)
  
  fragCounter <- 0
  
  for(i in 1:length(chrs)) {
  
    chr <- chrs[i]
  
    fragRanges <- ranges(matchPattern(patter=REmotif,subject=BScurrent[[chr]]))
    chrGR <- resize(GRanges(seqnames=chr,ranges=fragRanges),width=1,fix="end")
  
    N <- length(chrGR)
    start(chrGR)[1] <- 1
    end(chrGR)[N] <- sLengths[chr]+1
    start(chrGR)[2:N] <- end(chrGR)[1:(N-1)]+1
  
    chrGR$fragID <- (1:N)+fragCounter
    frags[[chr]] <- chrGR
  
    fragCounter <- fragCounter+N
  
  }
  
  if(writeToRDS) {
    
    saveRDS(frags,file=configOpt$fragsFile)
    
  } else {
    
    return(frags)
  }

}

createPartitionGR <- function(frags,partSize=8e6) {
  
  chrs <- names(frags)
  
  for(i in 1:length(chrs)) {
    
    chr <- chrs[i]
    
    S1 <- start(reduce(frags[[chr]]))
    E1 <- end(reduce(frags[[chr]]))
    
    S2 <- S1+floor(partSize/2)
    
    IRS1 <- seq(from=S1,to=E1-partSize,by=partSize)
    N <- length(IRS1)
    IRE1 <- c(IRS1[2:N]-1,E1)
    
    gR1 <- GRanges(chr,IRanges(IRS1,IRE1))
    
    IRS2 <- seq(from=S2,to=E1-partSize-floor(partSize/2),by=partSize)
    N <- length(IRS2)
    IRE2 <- c(IRS2[2:N]-1,E1-floor(partSize/2))
    gR2 <- GRanges(chr,IRanges(IRS2,IRE2))
    
    partGRchr <- sort(c(gR1,gR2))
    
    if(chr=="chr1") {
      
      partGR <- partGRchr
      
    } else {
      
      partGR <- suppressWarnings(c(partGR,partGRchr))
      
    }
    
  }
  
  partGR$partID <- paste0("part.",1:length(partGR))
  
  return(partGR)
  
}

getRecipPeaks <- function(peakHiCObj,anchorSize=10e3,vpSize=10e3,loopDF=NULL,loopFile=NULL) {
  #browser()
  vpsGR <- peakHiCObj$vpsGR
  
  if(!is.null(loopDF)){
    
    loops <- loopDF
    
  } else{
    
    if(is.null(loopFile)) {
      
      rdsFldr <- paste0(peakHiCObj$configOpt$projectFolder,"rds/")
      loopsFldr <- paste0(rdsFldr,"loops/")
      nReps <- sum(as.vector(peakHiCObj[["hic"]][["design"]][["HiCMap"]])==peakHiCObj$configOpt$hicCond)
      
      loopFile <- paste0(loopsFldr,peakHiCObj$name,"_GW_nReps_",nReps,"_peakHiC_wSize_",peakHiCObj$configOpt$peakCalls$wSize,"_qWr_",peakHiCObj$configOpt$peakCalls$qWr,"_alphaFDR_",peakHiCObj$configOpt$peakCalls$alphaFDR,"_loops.txt")
      
    }
    
    if(!file.exists(loopFile)) {
      
      message ("peakHiC loop file not found !")
      return(NULL)
      
    } else {
      
      loops <- fread(loopFile,stringsAsFactors=FALSE,sep="\t",header=TRUE)
      
    }
    
  }
  
  lx <- resize(GRanges(seqnames=loops$chr,IRanges(loops$vp_X1,loops$vp_X2)),width=anchorSize,fix="center")
  ly <- resize(GRanges(seqnames=loops$chr,IRanges(loops$maxV4CscorePos,loops$maxV4CscorePos)),width=anchorSize,fix="center")
  
  ovl <- findOverlaps(ly,resize(vpsGR,width=vpSize,fix="center"))
  
  loopPairs <- paste0(loops$id[ovl@from],".to.",as.vector(vpsGR$vpID)[ovl@to])
  revPairs <- paste0(as.vector(vpsGR$vpID)[ovl@to],".to.",loops$id[ovl@from])
  
  loops$loopID <- 1:nrow(loops)
  recipIDs <- unique(ovl@from[loopPairs%in%revPairs])
  
  loops$recip <- loops$loopID%in%recipIDs 
  
  recipLoops <- loops[recipIDs,]
  
  lx$id <- loops$id
  lx$loopID <- loops$loopID
  
  out <- list(loops=loops,lx=lx,ly=ly,recipLoops=recipLoops)
  
  return(out)
  
}


binGenome <- function(BSgenome,binSize=10e3){
  
  do.call( require, args=list( BSgenome ) )
  assign( 'genome', base::get(as.vector(BSgenome)))
  
  x<-seqinfo(genome)
  bins <- tileGenome(x, tilewidth=binSize, cut.last.tile.in.chrom=TRUE)
  
  return(bins)
  
}

getInsuBins <- function(peakHiCObj,binRes=2e3,scaleParam=200e3){
  
  BSgenome <- as.vector(strsplit(peakHiCObj$configOpt$BSgenome,split=" ")[[1]][[2]])
  
  do.call( require, args=list( BSgenome ) )
  assign( 'genome', base::get(as.vector(BSgenome)))
  
  x<-seqinfo(genome)
  
  chrs <- intersect(as.vector(unique(seqnames(peakHiCObj$partition$partGR))),names(seqlengths(genome)))
  
  S <- GRanges(seqnames=chrs,IRanges(rep(1,length(chrs)),rep(1,length(chrs))))
  E0 <- as.vector(seqlengths(genome)[chrs])
  E <- GRanges(seqnames=chrs,IRanges(E0,E0))
  
  edges <- c(S,E)
  
  insuBins <- tileGenome(x, tilewidth=binRes, cut.last.tile.in.chrom=TRUE)
  insuBins <- insuBins[as.vector(seqnames(insuBins))%in%chrs]
  
  D0 <- as.data.frame(distanceToNearest(insuBins,edges))
  D0 <- D0[D0$distance>scaleParam,]
  
  insuBins <- insuBins[D0$queryHits]
  insuBins$pos <- start(resize(insuBins,width=1,fix="center"))
  insuBins$binID <- paste0("bin.",1:length(insuBins))
  insuBins$partID <- NA
  
  partGR <- peakHiCObj$partition$partGR
  D1 <- distanceToNearest(resize(insuBins,width=1,fix="center"),resize(partGR,width=1,fix="center"))
  insuBins$partID[D1@from] <- partGR$partID[D1@to]
  
  return(insuBins)
  
}

getDI <- function(i,SX,SY) {
  
  xReads <- setdiff(SX[[i]],SY[[i]])
  yReads <- setdiff(SY[[i]],SX[[i]])
  
  B <- length(xReads)
  A <- length(yReads)
  
  E <- (A+B)/2
  
  di <- ((B-A)/abs(B-A))*((A-E)^2/E+(B-E)^2/E)
  
  return(di)
  
}

getINSU <- function(binID,SX,SY,PE1,PE2,bins,scaleParam=200e3) {
  
  IDX <- intersect(SX[[binID]],SY[[binID]])
  xPos <- start(PE1[IDX])
  yPos <- start(PE2[IDX])
  
  pos <- bins$pos[match(binID,bins$binID)]
  
  A <- length(IDX)
  B <- sum(yPos<pos)
  C <- sum(xPos>pos)
  
  return(log2((A-B-C)/A))
  
}

getPartitionDIs <- function(partID,peakHiCObj,trackID=NULL,binRes=2e3,binWidth=20e3,writeDIFile=FALSE) {
  
  if(is.null(trackID)) {
    
    trackID <- peakHiCObj$hic$design$trackID[1]
    
  }
  
  if(!is.null(peakHiCObj$data$insuBins)){
    
    insuBins <- readRDS(peakHiCObj$data$insuBins)
    
  } else {
    
    insuBins <- getInsuBins(peakHiCObj,binRes=binRes)
    
  }
  
  insuBins <- insuBins[insuBins$partID==partID]
  diBins <- resize(insuBins,width=binWidth,fix="center")
  
  PEdat <- getPartitionReads(partID=partID,trackID=trackID,peakHiCObj=peakHiCObj)
  
  olapX <- findOverlaps(diBins,PEdat$PE1)
  olapY <- findOverlaps(diBins,PEdat$PE2)
  SX <- split(olapX@to,olapX@from)
  SY <- split(olapY@to,olapY@from)
  
  binPresent <- intersect(names(SX),names(SY))
  
  SX <- SX[match(binPresent,names(SX))]
  SY <- SY[match(binPresent,names(SY))]
  
  dis <- sapply(1:length(binPresent),getDI,SX=SX,SY=SY)
  pos <- diBins$pos[as.numeric(binPresent)]
  chrs <- rep(as.vector(seqnames(diBins))[1],length(pos))
  out <- data.frame(chr=chrs,pos=pos,DI=dis,stringsAsFactors=FALSE)
  
  if(writeDIFile) {
    
    diFile <- paste0(peakHiCObj$configOpt$projectFolder,trackID,"_DI.txt")
    
    if(!file.exists(diFile)) {
      
      diHeader <- paste("chr","pos","DI",sep="\t")
      cat(c(diHeader,"\n"),file=diFile,sep="")
      
    }
    
    write.table(out,file=diFile,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
    
  } else {
    
    return(out)
    
  }
  
}

getBinPairs <- function(loops,bins,anchorSize=10e3) {
  
  lx <- resize(GRanges(seqnames=loops$chr,IRanges(loops$vp_X1,loops$vp_X2)),width=1,fix="center")
  ly <- resize(GRanges(seqnames=loops$chr,IRanges(loops$maxV4CscorePos,loops$maxV4CscorePos)),width=1,fix="center")
  
  ovl.lx <- findOverlaps(lx,bins)
  ovl.ly <- findOverlaps(ly,bins)
  
  loopIDs <- unique(intersect(ovl.lx@from,ovl.ly@from))
  loops <- loops[loopIDs,]
  
  ovl.lx <- ovl.lx[match(loopIDs,ovl.lx@from)]
  ovl.ly <- ovl.ly[match(loopIDs,ovl.ly@from)]
  
  NR.binID <- paste0(pmin(ovl.lx@to,ovl.ly@to),".",pmax(ovl.lx@to,ovl.ly@to))
  loops$NR.binID <- NR.binID
  
  return(loops)
  
}

getV4CRepReads <- function(repIDX,vpID,rds,hicCond) {
  
  return(rds[[vpID]][[hicCond]][[repIDX]]$GR$reads)
  
}

getV4CCovs <- function(vpID,rds,hicCond,nearCis=100e3) {
  
  nReps <- length(rds[[vpID]][[hicCond]])
  
  vpGR <- rds[[vpID]][[hicCond]][[1]][["vpGR"]]
  gR <- rds[[vpID]][[hicCond]][[1]]$GR
  
  if(! ( is.null(vpGR) | is.null(gR))) {
    
    nearCisOvlIDX <- findOverlaps(resize(vpGR,width=2*nearCis,fix="center"),gR)@to
    
    reads <- sapply(1:nReps,getV4CRepReads,vpID=vpID,rds=rds,hicCond=hicCond)
    
    if(!is.null(dim(reads))) { 
      
      totReads <- apply(reads,1,sum)
      
      if(length(nearCisOvlIDX)>1){
        
        nearCisCaptMean <- mean(apply(reads[nearCisOvlIDX,]>0,2,sum)/length(nearCisOvlIDX))
        nearCisCaptTot <- sum(totReads[nearCisOvlIDX]>0)/length(nearCisOvlIDX)
        nearCisCovMean <- mean(apply(reads[nearCisOvlIDX,],2,sum))
        nearCisCovTot <- sum(apply(reads[nearCisOvlIDX,],2,sum))
        
      } else {
        
        nearCisCaptMean <- 0
        nearCisCaptTot <- 0
        nearCisCovMean <- 0
        nearCisCovTot <- 0
      }
      
      plotCaptMean <- mean(apply(reads>0,2,mean))
      plotCaptTot <- sum(totReads>0)/length(totReads)
      
      plotCovMean <- mean(apply(reads,2,sum))
      plotCovTot <- sum(totReads)
      
      vpSize <- width(vpGR)
      
      return(c(nearCisCaptMean,nearCisCaptTot,nearCisCovMean,nearCisCovTot,plotCaptMean,plotCaptTot,plotCovMean,plotCovTot,vpSize))
      
    }
    
  } else {
    
    return(c(0,0,0,0,0,0,0,0,0))
    
  }
  
}

getPartV4CStats <- function(partID,peakHiCObj,hicCond=NULL,outFile=NULL) {
  
  if(is.null(outFile)) {
    
    outFile <- paste0(peakHiCObj$configOpt$projectFolder,"rds/",peakHiCObj$name,"_V4C_stats.txt")
    
  }
  
  if(is.null(hicCond)) {
    
    hicCond <- peakHiCObj$configOpt$hicCond
    
  }
  
  rdsFile <- paste0(peakHiCObj$configOpt$projectFolder,"rds/vpReads_",partID,".rds")
  
  if(file.exists(rdsFile)) {
    
    rds <- readRDS(file=rdsFile)
    vpIDs <- setdiff(names(rds),"configOpt")
    
    stats <- t(sapply(vpIDs,getV4CCovs,rds=rds,hicCond=hicCond))
    
  }
  
  write.table(stats,file=outFile,sep="\t",col.names=FALSE,row.names=TRUE,quote=FALSE,append=TRUE)
  
}

getV4CStatsByChr <- function(chr,peakHiCObj,hicCond=NULL,outFile=NULL) {
  
  if(is.null(outFile)) {
    
    outFile <- paste0(peakHiCObj$configOpt$projectFolder,"rds/",peakHiCObj$name,"_V4C_stats.txt")
    
  }
  
  if(is.null(hicCond)) {
    
    hicCond <- peakHiCObj$configOpt$hicCond
    
  }
  
  partGR <- subChr(gR=peakHiCObj$partition$partGR,chr=chr)
  
  rdsFolder <- paste0(peakHiCObj$configOpt$projectFolder,"rds/")
  partIDs <- gsub("vpReads_","",gsub("\\.rds","",dir(rdsFolder)[grep("vpReads_",dir(rdsFolder))]))
  
  partIDs <- intersect(partIDs,partGR$partID)
  
  if(length(partIDs)) {
    
    for(partID in partIDs){
      
      getPartV4CStats(partID=partID,peakHiCObj=peakHiCObj,hicCond=hicCond,outFile=outFile)
      
    }
    
  }
  
}

initExampleData <- function(baseFolder, pairixBinary) {

  fRDS <- paste0(baseFolder,"DATA/example_data/hg38_4DN_Rao_GM12878_peakHiC_example_peakHiCObj.rds")
  
  peakHiCObj <- readRDS(fRDS)
  peakHiCObj$configOpt$baseFolder <- baseFolder
  peakHiCObj$configOpt$pairixBinary <- pairixBinary
  peakHiCObj$configOpt$peakHiCObj <- paste0(baseFolder,"DATA/example_data/hg38_4DN_Rao_GM12878_peakHiC_example_peakHiCObj.rds")
  peakHiCObj$configOpt$projectFolder <- paste0(baseFolder,"RESULTS/Rao_4DN_GM12878_peakHiC_example/")
  peakHiCObj$configOpt$hicReadsFldr <- paste0(baseFolder,"DATA/example_data/")
  peakHiCObj$configOpt$fragsFile <- paste0(baseFolder,"DATA/hg38_MboI_fragsByChr.rds")
  saveRDS(peakHiCObj,file=peakHiCObj$configOpt$peakHiCObj)
  
  resultsFolder <- paste0(baseFolder,"RESULTS/")
  dataFolder <- paste0(resultsFolder,"Rao_4DN_GM12878_peakHiC_example/")
  
  if(!file.exists(resultsFolder)) {
  
    cmd <- paste0("mkdir ",resultsFolder)
    system(cmd)
    cmd <- paste0("mkdir ",dataFolder)
    system(cmd)
    
  }

}

peakHiCConf <- function(confFile){
  library(GenomicRanges)
  #browser()
  configOpt <- createConfig(confFile=confFile)
  designDF <- readHiCDesign(designFile=configOpt$hic$designFile)
  partGR <- readRDS(configOpt$genomePartition$partitionGR)
  
  peakHiCObj <- list()
  peakHiCObj$hic <- list()
  peakHiCObj$partition <- list()
  peakHiCObj$vpsGR <- NULL
  peakHiCObj$configOpt <- configOpt
  peakHiCObj$name <- configOpt$name
  peakHiCObj$hic$design <- designDF
  peakHiCObj$partition$partGR <- partGR
  
  vpData <- read.table(configOpt$VPsFile,stringsAsFactors=FALSE,header=TRUE)
  peakHiCObj$vpsGR <- createVPs(vpData=vpData,peakHiCObj=peakHiCObj)
  
  saveRDS(peakHiCObj,file=peakHiCObj$configOpt$peakHiCObj)
  
  message("\nPeakHiC configuration was read and results were stored in this RDS file :")
  message(peakHiCObj$configOpt$peakHiCObj,"\n")
  
}

exportLoops <- function(peakHiCObj,loopFile=NULL,makeNR=TRUE) {
  
  if(is.null(loopFile)) {
    
    rdsFldr <- paste0(peakHiCObj$configOpt$projectFolder,"rds/")
    loopsFldr <- paste0(rdsFldr,"loops/")
    wSize <- peakHiCObj$configOpt$peakCalls$wSize
    alphaFDR <- peakHiCObj$configOpt$peakCalls$alphaFDR
    qWr <- peakHiCObj$configOpt$peakCalls$qWr
    nReps <- nrow(peakHiCObj$hic$design)
    loopFile <- paste0(loopsFldr,peakHiCObj$name,"_GW_nReps_",nReps,"_peakHiC_wSize_",wSize,"_qWr_",qWr,"_alphaFDR_",alphaFDR,"_processed_loops.txt")
    
  }
  
  if(!file.exists(loopFile)) {
    
    stop("loopFile does not exist")
    
  }
  
  loopDF <- read.table(loopFile,header=TRUE,stringsAsFactors=FALSE)
  
  writeHICCUPs2D <- function(lx,ly,outFile,clr=NULL,width=NULL) {
    
    labs <- c("chr1","x1","x2","chr2","y1","y2","color")
    
    if(!is.null(width)) {
      
      lx <- resize(lx,width=width,fix="center")
      ly <- resize(ly,width=width,fix="center")
      
    } 
    
    out <- cbind(as.data.frame(lx)[,1:3],as.data.frame(ly)[,1:3])
    
    if(is.null(clr)) {
      
      out$color <- rep("0,255,255",nrow(out))
      
    } else {
      
      out$color <- rep(clr,nrow(out))
      
    }
    
    colnames(out) <- labs
    
    write.table(out,file=outFile,sep="\t",row.names=F,quote=F)
    
  }
  
  if(!is.null(loopDF$NR.binID) & makeNR) {
    
    loopDF <- loopDF[order(loopDF[[paste0(peakHiCObj$configOpt$hicCond,".covQ.norm")]],decreasing=TRUE),]
    loopDF <- loopDF[!duplicated(as.vector(loopDF$NR.binID)),]
    lxPos <- pmin(loopDF$vp_X1,loopDF$maxV4CscorePos)
    lyPos <- pmax(loopDF$vp_X1,loopDF$maxV4CscorePos)
    lx <- GRanges(loopDF$chr,IRanges(lxPos-5e3+1,lxPos+5e3))
    ly <- GRanges(loopDF$chr,IRanges(lyPos-5e3+1,lyPos+5e3))
    hiccupsFile <- paste0(loopsFldr,peakHiCObj$name,"_GW_nReps_",nReps,"_peakHiC_wSize_",wSize,"_qWr_",qWr,"_alphaFDR_",alphaFDR,"_loops_HICCUPS_format.txt")
    writeHICCUPs2D(lx=lx,ly=ly,outFile=hiccupsFile)
    
  }
  
}

combine.experiments <- function( data, num.exp = 0, vp.pos ){
  if(num.exp == 0){
    num.exp = data$num.exp
  }
  #create two element vector containing the viewpoint position
  #if only one viewpoint is given
  if(length(vp.pos) == 1){
    vp.pos <- c(vp.pos,vp.pos)
  }
  vp.pos <- sort(vp.pos)
  data.m <- data[[1]]
  #browser()
  for( i in 2:num.exp ){
    data.m <- base::merge(data.m, data[[i]], by=1)
  }
  #create the background model for the upstream regions
  data.bg <- data.m
  for( i in 1:num.exp ){
    data.bg[data.m[,1] < vp.pos[1],i+1] <- get.background(data.m[data.m[,1] < vp.pos[1],c(1,i+1)], vp.pos[1] )
  }
  #and for the downstream regions
  for( i in 1:num.exp ){
    data.bg[data.m[,1] > vp.pos[2],i+1] <- get.background(data.m[data.m[,1] > vp.pos[2],c(1,i+1)], vp.pos[2] )
  }
  #if two viewpoint fragments are given set the intervening fragments
  #to zero
  #set background to 1 to prevent NaN in the ratio
  if(vp.pos[1] != vp.pos[2]){
    for( i in 1:num.exp){
      data.m[data.m[,1] >= vp.pos[1] & data.m[,1] <= vp.pos[2],i+1] <- 0
      data.bg[data.m[,1] >= vp.pos[1] & data.m[,1] <= vp.pos[2],i] <- 1
    }
  }
  cbind(data.m, data.bg[,-1])
}

#perform pava regression and return the background regression line
get.background <- function( data, vp.pos, weight.factor=0, fractile=F){
  require(isotone)
  switched = FALSE
  weights <- (1:nrow(data))**weight.factor
  if(data[1,1] > vp.pos){
    data[,1] <- -data[,1] #reverse the sign to make the trend increasing
    switched = TRUE
    weights <- rev(weights)
  }
  #create the isotonic regression
  if(fractile){
    lm <- isotone::gpava(data[,1], data[,2], solver=weighted.fractile, weights=NULL, p=0.75)
  }else{
    lm <- isotone::gpava(data[,1], data[,2], solver=weighted.mean)
  }
  
  if(switched)
    pred.data <- data.frame( -lm$z, lm$x )
  else
    pred.data <- data.frame( lm$z, lm$x )
  
  pred.data[order(pred.data[,1]),2]
}

get.single.background <- function(data, num.exp = 1, vp.pos) {
  
  if (length(vp.pos) == 1) {
    vp.pos <- c(vp.pos, vp.pos)
  }
  vp.pos <- sort(vp.pos)
  
  data.bg <- data
  data.bg[data[, 1] < vp.pos[1], 2] <- get.background(data[data[, 1] < vp.pos[1], c(1, 2)], vp.pos[1])
  data.bg[data[, 1] > vp.pos[2], 2] <- get.background(data[data[, 1] > vp.pos[2], c(1, 2)], vp.pos[2])
  
  return(cbind(data, data.bg[, -1]))
  
}

makeRDSFolder <- function(rdsFldr){
  dir.create(rdsFldr, recursive = TRUE)
  dir.create(paste0(rdsFldr,"loops/"))
  dir.create(paste0(rdsFldr,"profiles/"))
  message("\nproject folder created at ", rdsFldr,"\n")
}


