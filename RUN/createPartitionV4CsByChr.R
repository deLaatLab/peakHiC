message('\n')
message( '#######################################################################################' )
message( '# peakHi-C - 2021 Geert Geeven and Valerio Bianchi - Hubrecht Institute ###############' )
message( '#######################################################################################' )
message('\n')
message('Creating virutal 4C profiles per partition per chromosome...')


#################################################################################################################
### PARSING THE INPUT ###########################################################################################
#################################################################################################################

if( !suppressMessages(require( "argparse", character.only = TRUE ) ) ) stop( "Package not found: argparse" )

options( scipen=999 )

# create parser object
parser <- ArgumentParser()

# specify our desired options
parser$add_argument('-chr', help='chromosome to analyze', metavar='ZZ,YY,AA,DD', required=TRUE )
parser$add_argument('-peakHiCObj', help='location of peakHiC peakHiCObject (rds file) with all parameters to run peakHiC', metavar='ZZ,YY,AA,DD', required=TRUE )

parser$add_argument('-wSize', type='integer', help='window size (# of fragments) for V4C profile running mean', default=NULL, metavar='XX' )
parser$add_argument('-vpSize', type='integer', help='# of fragments for V4C meta-viewpoint definition (pooling of contacts)', default=NULL, metavar='XX' )
parser$add_argument('-v4cSize', type='integer', help='genomic size (bp) of V4C contact profile', default=NULL, metavar='XX' )
parser$add_argument('-alphaFDR', type='double', help='false discovery rate (FDR) parameter', default=NULL, metavar='XX' )
parser$add_argument('-hicCond', type='double', help='name of the Hi-C data (cell type / condition) to analyze', default=NULL, metavar='XX' )
parser$add_argument('-baseFolder', type='character', help='path to the base (install) folder of peakHiC', default=NULL, metavar='XX' )
parser$add_argument('-projectFolder', type='character', help='path to the folder where peakHiC will store results of the analysis for this project', default=NULL, metavar='XX' )
parser$add_argument('-hicReadsFolder', type='character', help='path to the folder where the reads (Hi-C pairs files) are stored', default=NULL, metavar='XX' )

args <- parser$parse_args()

chr <- args$chr
peakHiCObjFile <- args$peakHiCObj

wSize <- args$wSize
vpSize <- args$vpSize
alphaFDR <- args$alphaFDR
v4cSize <- args$v4cSize
hicCond <- args$hicCond
baseFolder <- args$baseFolder
projectFolder <- args$projectFolder
hicReadsFolder <- args$hicReadsFolder


#################################################################################################################
### Load packages ###############################################################################################
#################################################################################################################

message( paste0( '> loading R libraries and functions..' ) )
if( !suppressMessages(require( "GenomicRanges", character.only=TRUE ) ) ) stop( "Package not found: GenomicRanges" )
if( !suppressMessages(require( "zoo", character.only=TRUE ) ) ) stop( "Package not found: zoo" )
if( !suppressMessages(require( "parallel", character.only=TRUE ) ) ) stop( "Package not found: parallel" )
if( !suppressMessages(require( "data.table", character.only=TRUE ) ) ) stop( "Package not found: data.table" )

#################################################################################################################
### Load object with global peakHiC settings ####################################################################
#################################################################################################################
message( paste0( '> loading peakHiCObj..' ) )

#chr <- 'chr5'
#peakHiCObjFile <-"~/data/Leducq/peakHiCObj/hg38_Leducq_LA_LV_12kbbins_chr5.rds" 
chr <- args$chr 
peakHiCObjFile <- args$peakHiCObj 
peakHiCObj <- tryCatch({readRDS(peakHiCObjFile)}, error=function(cond){message('\npeakHiC ERROR: path to peakHiCObj incorrect\n')})

#override input from peakHiCObj if given as command line input
if(!is.null(wSize)) {peakHiCObj$configOpt$V4C$wSize <- wSize}
if(!is.null(wSize)) {peakHiCObj$configOpt$peakCalls$wSize <- wSize}
if(!is.null(vpSize)) {peakHiCObj$configOpt$V4C$vpSize <- vpSize}
if(!is.null(v4cSize)) {peakHiCObj$configOpt$V4C$v4cSize <- v4cSize}
if(!is.null(alphaFDR)) {peakHiCObj$configOpt$peakCalls$alphaFDR <- alphaFDR}
if(!is.null(hicCond)) {peakHiCObj$configOpt$hicCond <- hicCond}
if(!is.null(baseFolder)) {peakHiCObj$configOpt$baseFolder <- baseFolder}
if(!is.null(projectFolder)) {peakHiCObj$configOpt$projectFolder <- projectFolder}
if(!is.null(hicReadsFolder)) {peakHiCObj$configOpt$hicReadsFldr <- hicReadsFolder}

#################################################################################################################
### Load R functions and files ##################################################################################
#################################################################################################################
#source the file containing functions
sourceFile <- paste0(peakHiCObj$configOpt$baseFolder,"R/peakHiC_functions.R")
tryCatch({source(sourceFile)}, error=function(cond){message("\npeakHiC ERROR: baseFolder is incorrect, cannot call peakHiC_functions.R\n")} )

#load frag file
frags <- tryCatch({readRDS(peakHiCObj$configOpt$fragsFile)[[chr]]}, error=function(cond){message('\npeakHiC ERROR: fragsFile is indicated incorrectly\n')})

#create output folder
rdsFldr <- paste0(peakHiCObj$configOpt$projectFolder,"rds/")
if(!file.exists(rdsFldr)) { makeRDSFolder(rdsFldr) }

#main function to be called in parallel 
getReads.PartID <- function(partID) {
  
  source(sourceFile)
	fRDS <- paste0(rdsFldr,"vpReads_",partID,".rds")
	message( paste0( '> analyzing partition: ' , partID) )	

	tryCatch({
		vpReads <- getPeakHiCData(partID=partID,frags=frags,peakHiCObj=peakHiCObj)
		saveRDS(vpReads,file=fRDS)
			}, error=function(cond) { message(paste0("generation of V4Cs failed for partID ", partID))}) #conditionMessage()
   
}


###################################################################################
###Parallel call###################################################################
###################################################################################
#set n.o. cores
inputCores <- peakHiCObj$configOpt$nCores
detectedCores <- detectCores()
n.cores <- ifelse(inputCores < detectedCores, inputCores, detectedCores)
if(!is.integer(n.cores)){message("\nincorrect format for number of cores\n")}
message("\nRunning peakHiC on ", n.cores, " cores.\n")

cl = makeCluster(n.cores, outfile="")
clusterEvalQ(cl, c(suppressPackageStartupMessages({ 
  library("GenomicRanges") 
  library("data.table") 
  library("zoo") })) )

clusterExport(cl=cl, varlist=c("peakHiCObj", "sourceFile", "frags", "rdsFldr"), envir=environment())
#clusterEvalQ(cl, sessionInfo())


#partitions to run the call over
partIDs <- unique(subChr(peakHiCObj[["vpsGR"]],chr)$partID)
if(is.null(partIDs)){print("ERROR, incorrectly incoorporated viewpoints into peakHiCObj")}

##### call the function in parallel for each partID and write in rdsFldr'
#lapply(partIDs[1], FUN = getReads.PartID)
parLapply(cl, partIDs, fun = getReads.PartID)

stopCluster(cl)
message('>>>> DONE <<<<')
q("no")

