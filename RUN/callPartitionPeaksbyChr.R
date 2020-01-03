message('\n')
message( '#######################################################################################' )
message( '# peakHi-C - Copyright (C) 2018 Geert Geeven and Valerio Bianchi - Hubrecht Institute #' )
message( '#######################################################################################' )
message('\n')

#################################################################################################################
### DEFINE INPUT ARGS ###########################################################################################
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
parser$add_argument('-minDist', type='integer', help='minimum distance between vp and anchor in peak calling', default=NULL, metavar='XX' )
parser$add_argument('-alphaFDR', type='double', help='false discovery rate (FDR) parameter', default=NULL, metavar='XX' )
parser$add_argument('-qWr', type='double', help='peakHiC qWr ratio statistic parameter', default=NULL, metavar='XX' )
parser$add_argument('-hicCond', type='double', help='name of the Hi-C data (cell type / condition) to analyze', default=NULL, metavar='XX' )
parser$add_argument('-baseFolder', type='character', help='path to the base (install) folder of peakHiC', default=NULL, metavar='XX' )
parser$add_argument('-projectFolder', type='character', help='path to the folder where peakHiC will store results of the analysis for this project', default=NULL, metavar='XX' )
parser$add_argument('-hicReadsFolder', type='character', help='path to the folder where the reads (Hi-C pairs files) are stored', default=NULL, metavar='XX' )

#################################################################################################################
### Parse cmd line args #########################################################################################
#################################################################################################################

args <- parser$parse_args()

chr <- args$chr
peakHiCObjFile <- args$peakHiCObj

wSize <- args$wSize
vpSize <- args$vpSize
alphaFDR <- args$alphaFDR
minDist <- args$minDist
qWr <- args$qWr
v4cSize <- args$v4cSize
hicCond <- args$hicCond
baseFolder <- args$baseFolder
projectFolder <- args$projectFolder
hicReadsFolder <- args$hicReadsFolder

#################################################################################################################
### Load object with global peakHiC settings ####################################################################
#################################################################################################################

message( paste0( '> loading peakHiCObj..' ) )

peakHiCObj <- readRDS(peakHiCObjFile)

#################################################################################################################
### override global settings with specified cmd line args #######################################################
#################################################################################################################

if(!is.null(wSize)) {peakHiCObj$configOpt$V4C$wSize <- wSize}
if(!is.null(wSize)) {peakHiCObj$configOpt$peakCalls$wSize <- wSize}
if(!is.null(vpSize)) {peakHiCObj$configOpt$V4C$vpSize <- vpSize}
if(!is.null(v4cSize)) {peakHiCObj$configOpt$V4C$v4cSize <- v4cSize}
if(!is.null(minDist)) {peakHiCObj$configOpt$peakCalls$minDist <- minDist}
if(!is.null(alphaFDR)) {peakHiCObj$configOpt$peakCalls$alphaFDR <- alphaFDR}
if(!is.null(qWr)) {peakHiCObj$configOpt$peakCalls$qWr <- qWr}
if(!is.null(hicCond)) {peakHiCObj$configOpt$hicCond <- hicCond}
if(!is.null(baseFolder)) {peakHiCObj$configOpt$baseFolder <- baseFolder}
if(!is.null(projectFolder)) {peakHiCObj$configOpt$projectFolder <- projectFolder}
if(!is.null(hicReadsFolder)) {peakHiCObj$configOpt$hicReadsFldr <- hicReadsFolder}

#################################################################################################################
### Load R functions ############################################################################################
#################################################################################################################

sourceFile <- paste0(peakHiCObj$configOpt$baseFolder,"R/peakHiC_functions.R")

#################################################################################################################
### Load packages and R functions ###############################################################################
#################################################################################################################

message( paste0( '> loading R libraries and functions..' ) )

if( !suppressMessages(require( "GenomicRanges", character.only=TRUE ) ) ) stop( "Package not found: GenomicRanges" )
if( !suppressMessages(require( "zoo", character.only=TRUE ) ) ) stop( "Package not found: zoo" )
if( !suppressMessages(require( "peakC", character.only=TRUE ) ) ) stop( "Package not found: peakC" )
if( !suppressMessages(require( "data.table", character.only=TRUE ) ) ) stop( "Package not found: data.table" )

source(sourceFile)

#################################################################################################################
### Run peakHiC script ##########################################################################################
#################################################################################################################

frags <- readRDS(peakHiCObj$configOpt$fragsFile)[[chr]]

vpsGR <- peakHiCObj[["vpsGR"]]
ids <- unique(subChr(vpsGR,chr)$partID)

rdsFldr <- paste0(peakHiCObj$configOpt$projectFolder,"rds/")

for(i in 1:length(ids)) {

	partID <- ids[i]
	fRDS <- paste0(rdsFldr,"vpReads_",partID,".rds")

	if(file.exists(fRDS)){
	
		message( paste0( '> calling peaks with peakHiC in partition ..' , partID) )	

		tryCatch({
	
			getPartitionPeaks(partID=partID,peakHiCObj=peakHiCObj,writePeaksFile=TRUE)

		}, error=function(e) { message(paste0("peak calling failed for partID ", ids[i])) })
   }
}

message( paste0( '>>>> DONE <<<<' ) )

q("no")
