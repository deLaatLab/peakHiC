message('\n')
message( '#######################################################################################' )
message( '# peakHi-C - Copyright (C) 2019 Geert Geeven and Valerio Bianchi - Hubrecht Institute #' )
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

parser$add_argument('-loopsFile', type='character', help='input file containing peakHiC loops to be analyzed', default=NULL, metavar='XX' )
parser$add_argument('-outFile', type='character', help='window size (# of fragments) for V4C profile running mean', default=NULL, metavar='XX' )
parser$add_argument('-baseFolder', type='character', help='path to the base (install) folder of peakHiC', default=NULL, metavar='XX' )
parser$add_argument('-anchorSize', type='integer', help='size of loop anchors for coverage calculation', default=10000, metavar='XX' )

#################################################################################################################
### Parse cmd line args #########################################################################################
#################################################################################################################

args <- parser$parse_args()

chr <- args$chr
peakHiCObjFile <- args$peakHiCObj

loopsFile <- args$loopsFile
outFile <- args$outFile
baseFolder <- args$baseFolder
anchorSize <- args$anchorSize

#################################################################################################################
### Load object with global peakHiC settings ####################################################################
#################################################################################################################

if( !file.exists(peakHiCObjFile) ) stop( "file with peakHiCObj not found !" )

message( paste0( '> loading peakHiCObj..' ) )

peakHiCObj <- readRDS(peakHiCObjFile)

#################################################################################################################
### override global settings with specified cmd line args #######################################################
#################################################################################################################

hicCond <- peakHiCObj$configOpt$hicCond

designMat <- peakHiCObj[["hic"]][["design"]]
hicTracksByCondition <- split(as.vector(designMat$trackID),designMat$condID)
tracks <- hicTracksByCondition[[hicCond]]

nReps <- length(tracks)

rdsFldr <- paste0(peakHiCObj$configOpt$projectFolder,"rds/")
loopsFldr <- paste0(rdsFldr,"loops/")
	
if(is.null(loopsFile)) {
		
	loopsFile <- paste0(loopsFldr,peakHiCObj$name,"_GW_nReps_",nReps,"_peakHiC_wSize_",peakHiCObj$configOpt$peakCalls$wSize,"_qWr_",peakHiCObj$configOpt$peakCalls$qWr,"_alphaFDR_",peakHiCObj$configOpt$peakCalls$alphaFDR,"_loops.txt")
	
}

if( !file.exists(loopsFile) ) stop( "loopsFile not found !" )

if(!is.null(baseFolder)) {peakHiCObj$configOpt$baseFolder <- baseFolder}

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

loops <- fread(file=loopsFile,sep="\t",header=TRUE,stringsAsFactors=FALSE)

if(is.null(loops$loopID)) {loops$loopID <- 1:nrow(loops)}

getLoopCovbyChr(chr=chr,loops=loops,peakHiCObj=peakHiCObj,anchorSize=anchorSize,loopCovFile=outFile,writePeaksFile=TRUE) 

message( paste0( '>>>> DONE <<<<' ) )

q("no")
