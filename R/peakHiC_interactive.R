createVPs <- function(vpData, peakHiCObj, fragsFile=NULL) {
  
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
    distMap <- distanceToNearest(vpsGR,resize(partGR,width=1,fix="center"))
    
    vpsGR$partID <- "no.part"
    vpsGR$partID[distMap@from] <- partGR$partID[distMap@to]
    vpsGR <- sort(vpsGR[vpsGR$partID!="no.part"])
    
    fragsList <- readRDS(fragsFile)
    chrs <- intersect(names(fragsList),as.vector(seqnames(vpsGR)))
    
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

getVPs <- function(partID,peakHiCObj) {
  
  if(is.null(peakHiCObj[["vpsGR"]])) {
    
    out <- NULL
    
  } else {
    
    vpsGR <- peakHiCObj[["vpsGR"]]
    out <- vpsGR[vpsGR$partID==partID]
    
  }
  
  return(out)
}

exportV4Cbw <- function(vpID,peakHiCObj,trackFldr=NULL) {

  if( !suppressMessages(require( "rtracklayer", character.only=TRUE ) ) ) stop( "Bioconductor package not found: rtracklayer" )
  if(is.null(trackFldr)){
    
    if(is.null(peakHiCObj$configOpt$trackFldr)) {
      
      trackFldr <- paste0(peakHiCObj$configOpt$projectFolder,"TRACKS/")
      
    } else {
      
      trackFldr <- peakHiCObj$configOpt$trackFldr
      
    }
  
  }
  
  if(!file.exists(trackFldr)){
    dir.create(trackFldr)
  }
  
  vpIDX <- match(vpID,peakHiCObj$vpsGR$vpID)
  
  if(is.na(vpIDX)) {
    stop("vpID not found in vpsGR object")
  }
  
  nextVP <- peakHiCObj$vpsGR[vpIDX]
  fRDS <- paste0(peakHiCObj$configOpt$projectFolder,"rds/profiles/",nextVP$partID,".rds")
  
  if(!file.exists(fRDS)){
    stop("V4C data not found for partition containing requested vpID, please run createPartitionV4CsByChr.R first")  
  }
  
  hicCond <- peakHiCObj$configOpt$hicCond
  V4CData <- readRDS(fRDS)
  vpV4C <- V4CData[[hicCond]][[vpID]]
  
  BSname <- strsplit(peakHiCObj$configOpt$BSgenome,split=" ")[[1]][2]
  do.call(require,args=list(BSname))
  assign( 'BScurrent', base::get( BSname ) )
  
  gR <- GRanges(seqnames=seqnames(vpV4C),ranges=ranges(vpV4C),seqlengths=seqlengths(BScurrent),seqinfo=seqinfo(BScurrent))
  gR$score <- vpV4C$normV4C
  
  bwFile <- paste0(trackFldr,peakHiCObj$configOpt$hicCond,"_",vpID,"_V4C_track.bw")
  export(object=gR,con=bwFile)
  
}

getV4CData <- function(vpID,peakHiCObj,configOpt=NULL,wSize=21,alphaFDR=0.1,qWr=2.0,minDist=20e3,storeVPReads=FALSE) {
  
  if(is.null(configOpt)) {
    
    configOpt <- peakHiCObj$configOpt
    
  }
  
  rdsFldr <- paste0(configOpt$projectFolder,"rds/")
  v4cFldr <- paste0(configOpt$projectFolder,"rds/profiles/")
  
  fragsFile <- configOpt$fragsFile
  hicCond <- configOpt$hicCond
  
  vpsGR <- peakHiCObj[["vpsGR"]]
  vpIDX <- match(vpID,vpsGR$vpID)
  
  if(is.na(vpIDX)) {
    
    stop(paste0("Viewpoint with vpID: ",vpID," not found."))
    
  }
  
  partID <- vpsGR$partID[vpIDX]
  
  vpPos <- start(vpsGR[vpIDX])
  vpChr <- as.vector(seqnames(vpsGR[vpIDX]))[1]
  
  v4cRDS <- paste0(v4cFldr,partID,".rds")

  if(file.exists(v4cRDS)) {
    
    partData <- readRDS(v4cRDS)
    vpData <- partData[[hicCond]][[vpID]]
    
    if(!is.null(vpData)) {
      
      v4cDat <- data.frame(pos=start(vpData),normV4C=vpData$normV4C,chr=as.vector(seqnames(vpData)))
      return(v4cDat)
      
    } else {
      
      stop(paste0("V4C profile data not found for vp: ",vpID,"\nPlease run peakHiC pipeline first to create V4C data."))
      
    }
    
  } else {
    
    stop(paste0("V4C profile data not found for vp: ",vpID,"\nPlease run peakHiC pipeline first to create V4C data."))
    
  }
  
}

v4cPlot <- function(vpID,peakHiCObj,showLoops=FALSE,overlapGRs=NULL,loopFile=NULL,loopYlim=c(0.4,0.8),xdiv=1e6,...){

  pDat <- getV4CData(vpID=vpID,peakHiCObj=peakHiCObj)
  plot(x=pDat$pos/xdiv,y=pDat$normV4C,type="h",frame.plot=FALSE,ylab="V4C",xlab="pos (Mb)",...)
  
  if(showLoops){
      
      ylimCurrent <- par("usr")[3:4]
      loops <- getOvLoops(peakHiCObj=peakHiCObj,overlapGRs=overlapGRs,loopFile=loopFile)
      
      addLoops(peakHiCObj = peakHiCObj, loops=loops, ylimCurrent=ylimCurrent,loopYlim=loopYlim)
    
  }

}

addLoops <- function(peakHiCObj,loops,anchorPlotSize=20e3,xdiv=1e6,clr="darkgray",ylimCurrent=NULL,loopYlim=c(0.4,0.8)){
  
  if(is.null(ylimCurrent)) {
    
    ylimCurrent <- par("usr")[3:4]
  
  }
  
  lx <- resize(loops$lx,width=anchorPlotSize,fix="center")
  ly <- resize(loops$ly,width=anchorPlotSize,fix="center")
  
  if(length(lx)>0){
    
    totalRange <- c((ylimCurrent[2]-ylimCurrent[1])*loopYlim[1],(ylimCurrent[2]-ylimCurrent[1])*loopYlim[2])
    loopYms <- seq(from=totalRange[1],to=totalRange[2],length.out=length(lx))
    yOffSet <- (loopYms[2]-loopYms[1])/3
    
    for( x in 1:length( lx ) ) {
      
      plotRegions( lx[x], yrange=c( loopYms[x]-yOffSet, loopYms[x]+yOffSet ), col=clr )
      plotRegions( ly[x], yrange=c( loopYms[x]-yOffSet, loopYms[x]+yOffSet ), col=clr )
      
      minX <- ( min( start( lx[x] ), start( lx[x] ), end( ly[x] ), end( ly[x] ) )+5e3 ) / 1e6
      maxX <- ( max( start( lx[x] ), start( lx[x] ), end( ly[x] ), end( ly[x] ) )-5e3 ) / 1e6
      
      segments( x0=minX, x1=maxX, y0=loopYms[x], lty=2 )
      
    }
    
  }

}

plotRegions <- function( ranges, yrange, xdiv=1e6, plotRanges=NULL, col="black" ) {
  
  if(is.null(plotRanges)) {
    
  } else {
    
    ovlIdx <- findOverlaps(ranges,plotRanges)@queryHits
    ranges <- ranges[ovlIdx]
    
  }
  
  if(length(ranges)){
    
    rect( xleft=start(ranges)/xdiv, ybottom=yrange[1], xright=end(ranges)/xdiv, ytop=yrange[2], col=col, border=FALSE )
    
  }
  
}

getOvLoops <- function(peakHiCObj,overlapGRs=NULL,loopFile=NULL) {
  
  require(GenomicRanges)
  
  out <- NULL
  
  if(is.null(loopFile)) {
    
    rdsFldr <- paste0(peakHiCObj$configOpt$projectFolder,"rds/")
    loopsFldr <- paste0(rdsFldr,"loops/")
    wSize <- peakHiCObj$configOpt$peakCalls$wSize
    alphaFDR <- peakHiCObj$configOpt$peakCalls$alphaFDR
    qWr <- peakHiCObj$configOpt$peakCalls$qWr
    nReps <- nrow(peakHiCObj$hic$design)
    loopFile <- paste0(loopsFldr,peakHiCObj$name,"_GW_nReps_",nReps,"_peakHiC_wSize_",wSize,"_qWr_",qWr,"_alphaFDR_",alphaFDR,"_processed_loops.txt")
  
  }
  
  if(file.exists(loopFile)){
    
    loopDF <- read.table(file=loopFile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
    rankVar <- paste0(peakHiCObj$configOpt$hicCond,".covQ.norm")
    
    if(!is.null(loopDF[[rankVar]])){
      
      loopDF <- loopDF[order(loopDF[[rankVar]],decreasing=TRUE),]
      loopDF <- loopDF[!duplicated(loopDF$NR.binID),]
      
    }
    
    lx <- resize(GRanges(seqnames=loopDF$chr,IRanges(loopDF$vp_X1,loopDF$vp_X2)),width=10e3,fix="center")
    ly <- resize(GRanges(seqnames=loopDF$chr,IRanges(loopDF$maxV4CscorePos,loopDF$maxV4CscorePos)),width=10e3,fix="center")
    
    IDX <- 1:length(lx)
    
    if(!is.null(overlapGRs)) {  
    
      IDX <- which(countOverlaps(lx,overlapGRs)>0|countOverlaps(ly,overlapGRs)>0)
    
    }
    
    if(length(IDX)>0){
    
      out <- list(lx=resize(lx[IDX],width=10e3,fix="center"),ly=resize(ly[IDX],width=10e3,fix="center"))
    
    }
  
  } else {
      
      stop("Loopfile not found. Please check the path to the peakHiC loop file.")
  }
  
  return(out)
  
}

getLoopFile <- function(peakHiCObj) {
  
  rdsFldr <- paste0(peakHiCObj$configOpt$projectFolder,"rds/")
  loopsFldr <- paste0(rdsFldr,"loops/")
  wSize <- peakHiCObj$configOpt$peakCalls$wSize
  alphaFDR <- peakHiCObj$configOpt$peakCalls$alphaFDR
  qWr <- peakHiCObj$configOpt$peakCalls$qWr
  nReps <- nrow(peakHiCObj$hic$design)
  loopFile <- paste0(loopsFldr,peakHiCObj$name,"_GW_nReps_",nReps,"_peakHiC_wSize_",wSize,"_qWr_",qWr,"_alphaFDR_",alphaFDR,"_processed_loops.txt")
  
  return(loopFile)

}

exportLoops <- function(peakHiCObj,loopFile=NULL,outFilePrefix=NULL,makeNR=TRUE,vpFilter=NULL,loopType=NULL) {
  
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
    loopDF[order(loopDF$chr,loopDF$vp_X1,loopDF$anchor_X1),]
    
    loopTypes <- peakHiCObj$vpsGR$type[match(loopDF$id,peakHiCObj$vpsGR$vpID)]
    lxPos <- pmin(loopDF$vp_X1,loopDF$maxV4CscorePos)
    lyPos <- pmax(loopDF$vp_X1,loopDF$maxV4CscorePos)
    lx <- GRanges(loopDF$chr,IRanges(lxPos-5e3+1,lxPos+5e3))
    ly <- GRanges(loopDF$chr,IRanges(lyPos-5e3+1,lyPos+5e3))
    
    if(!is.null(vpFilter)) {
      
      ovlLoops <- countOverlaps(lx,vpFilter)>0&countOverlaps(ly,vpFilter)>0
      loopDF <- loopDF[ovlLoops,]

    }
    
    # filter loops by loopType (e.g. CTCF/TSS/H3K27ac)
    
    if(!is.null(loopType)) {
      
      loopDF <- loopDF[loopTypes%in%loopType,]

    }
    
    loopTypes <- peakHiCObj$vpsGR$type[match(loopDF$id,peakHiCObj$vpsGR$vpID)]
    lxPos <- pmin(loopDF$vp_X1,loopDF$maxV4CscorePos)
    lyPos <- pmax(loopDF$vp_X1,loopDF$maxV4CscorePos)
    lx <- GRanges(loopDF$chr,IRanges(lxPos-5e3+1,lxPos+5e3))
    ly <- GRanges(loopDF$chr,IRanges(lyPos-5e3+1,lyPos+5e3))
    
    if(is.null(outFilePrefix)) {
      
      outFilePrefix <- paste0(loopsFldr,peakHiCObj$name,"_GW_nReps_",nReps,"_peakHiC_wSize_",wSize,"_qWr_",qWr,"_alphaFDR_",alphaFDR)
      
    }
    
    hiccupsFile <- paste0(outFilePrefix,"_loops_HICCUPS_format.txt")
    writeHICCUPs2D(lx=lx,ly=ly,outFile=hiccupsFile)
    
    UCSCFile <- paste0(outFilePrefix,"_loops_UCSC_interact.bed")
    
    trackPrefix <- "track type=interact name=\""
    trackOptions <- " color=255,153,0 maxHeightPixels=200:100:50 visibility=full "
    browserCmd <- paste0("browser position chr1:42,127,739-44,067,118")
    headerLine <- paste0(trackPrefix,peakHiCObj$name,"\" description=\"peakHiC loops\"",trackOptions)
    
    cat(browserCmd,file=UCSCFile,sep="\n")
    cat(headerLine,file=UCSCFile,sep="\n",append=TRUE)
    
    UCSC_Interact_DF <- data.frame(chrom = loopDF$chr,
                                              allStart = lxPos,
                                              allEnd = lyPos,
                                              name = loopDF$id,
                                              score = rep("0",nrow(loopDF)),
                                              value = rep("1",nrow(loopDF)),
                                              exp = loopTypes,
                                              color = rep("255,153,0",nrow(loopDF)),
                                              chromLower = loopDF$chr,
                                              startLower = lxPos-1,
                                              endLower = lxPos,
                                              name2 = loopDF$id,
                                              strand = rep(".",nrow(loopDF)),
                                              chromUpper = loopDF$chr,
                                              startUpper = lyPos-1,
                                              endUpper = lyPos,
                                              name3 = loopDF$id,
                                              strand2 = rep(".",nrow(loopDF)),stringsAsFactors=FALSE)
    
    write.table(UCSC_Interact_DF,file=UCSCFile,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE,quote=FALSE)
    
  }
  
}
