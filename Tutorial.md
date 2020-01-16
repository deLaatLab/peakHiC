# peakHiC - interactive data analysis and export

If you have followed the steps in the README.md Markdown document and have installed the peakHiC pipeline and have run the commands to create V4C profiles and loops for the example data, the resulting files should now be available for further analysis. Here we will briefly describe how to visualize V4C plots from viewpoints of interest in R and how to export loops and V4C BigWig tracks for further analysis using e.g. the UCSC (https://genome.ucsc.edu/) or IGV (https://igv.org/) genome browsers or the Juicebox (https://github.com/aidenlab/Juicebox/wiki) tool developed by the lab of Erez Lieberman Aiden. 

The first step is to open the R console and load the peakHiC interactive functions from this repository. Again make sure that you set the __baseFolder__ variable in the code below correctly to match your local installation.

```{r source}
baseFolder <- "/home/geert/localdev/github/peakHiC/"
sourceFile <- paste0(baseFolder,"R/peakHiC_interactive.R")
source(sourceFile)
```

Anytime we want to analyze data with peakHiC, we also need to load a peakHiCObj that contains all configuration data for a specific peakHiC dataset and analysis run. Here, we will use the peakHiC example data, analyzing Hi-C from GM12878 cells limited to a ~4Mb section (181 viewpoints) on chr1 (hg38). 

```{r source}
objFile <- paste0(baseFolder,"DATA/example_data/hg38_4DN_Rao_GM12878_peakHiC_example_peakHiCObj.rds")
peakHiCObj <- readRDS(objFile)
```

Viewpoints can be defined from any locus of interest in the genome. The example data contains 181 viewpoints from gene promoters (TSS), CTCF bound site in GM12878 (ENCODE ChIP-Seq data) and sites with enrichment for H3K27ac in GM12878. For the example data, these viewpoints are already stored in a GenomicRanges object. You can view them by typing

```{r source}
peakHiCObj$vpsGR
```

in the R console. To add your own viewpoints, you need to create a txt file which defines them. 

* Viewpoint file
  * peakHiC viewpoints are organized in a viewpoint file. Columns in this file specify the genomic coordinates of each viewpoint (row in the file). Additionally we need to specify a unique vpID, a name and a type (e.g. TSS, CTCF site) : 

| chr   | vpPos     | vpID              | name              | type |
|-------|-----------|-------------------|-------------------|------|
| chr1  | 42846618  | CTCF_ENCODE_14013 | CTCF_ENCODE_14013 | CTCF |
| chr1  | 42931204  | CTCF_ENCODE_25928 | CTCF_ENCODE_25928 | CTCF |

**Table 1.** Example of a peakHiC viewpoint file which defines genomic loci from which V4C profiles will be created.

The R code below will read this file and replace the vpsGR object in peakHiCObj, so that these viewpoints can be analyzed. Based on the genomic partition and restriction fragments defined in the peakHiC object, the __createVPs__ function will assign a partID and fragID to each viewpoint. 

```{r source}
exampleVPs <- paste0(baseFolder,"DATA/example_data/hg38_4DN_Rao_GM12878_peakHiC_example_VPs.txt")
vpData <- read.table(exampleVPs,header=TRUE,stringsAsFactors=FALSE)
peakHiCObj$vpsGR <- createVPs(vpData=vpData,peakHiCObj=peakHiCObj)
```

If for instance we want to plot the V4C profile for one of the peakHiC viewpoints from **Table 1**, we can use the following R code 

```{r source}
vpID <- "CTCF_ENCODE_25928"
v4cPlot(vpID=vpID,peakHiCObj=peakHiCObj,ylim=c(0,2))
```

This will prompt R to open a display with the resulting plot, which should look like this:

![peakHiC BigWig track in IGV](https://github.com/deLaatLab/peakHiC/raw/master/tutorial/peakHiC_example_v4cPlot_R_CTCF_ENCODE_25928.png)

**Figure 1.** Visualization of peakHiC V4C profile interactively in an R session

To add loops to this plot, we can load the entire set of peakHiC loops, or make a selection of loops only relevant for this locus. You can specify loops overlapping any genomic regions of interest (such as promoters / enhancers / CTCF sites), but here we will just plot ALL loops overlapping with 2 example viewpoints:

```{r source}
ids <- c("CTCF_ENCODE_25928","CTCF_ENCODE_35219")
overlapGRs <- peakHiCObj$vpsGR[match(ids,peakHiCObj$vpsGR$vpID)]
par(mfrow=c(2,1))
clr <- "#ff9900"
v4cPlot(vpID=ids[1],peakHiCObj=peakHiCObj,showLoops=TRUE,overlapGRs=overlapGRs,ylim=c(0,2),col=clr)
xlim=par("usr")[1:2]
v4cPlot(vpID=ids[2],peakHiCObj=peakHiCObj,showLoops=TRUE,overlapGRs=overlapGRs,xlim=xlim,ylim=c(0,2),col=clr)
```

![peakHiC BigWig track in IGV](https://github.com/deLaatLab/peakHiC/raw/master/tutorial/peakHiC_example_v4cPlot_R_CTCF_VPs_with_loops.png)

**Figure 2.** Visualization of V4C profiles of 2 example CTCF viewpoints together with peakHiC loops

We can also export the peakHiC V4C tracks as BigWig tracks, which can be added to a trackHub for visualization in the UCSC browser or we can directly load them into IGV. You first need to specify a folder to write the tracks to. Make sure this folder exists and you have permission to write files to it. The code below will export 2 example tracks to this folder.

```{r source}
trackFldr <- "/home/geert/localdev/github/peakHiC/RESULTS/Rao_4DN_GM12878_peakHiC_example/TRACKS/"
exportV4Cbw(vpID="CTCF_ENCODE_14013",peakHiCObj=peakHiCObj,trackFldr=trackFldr)
exportV4Cbw(vpID="CTCF_ENCODE_25928",peakHiCObj=peakHiCObj,trackFldr=trackFldr)
```

These tracks are ready to be imported into UCSC or the IGV genome browser. Below is a screenshot where we loaded the exported BigWig tracks into the IGV browser.

![peakHiC BigWig track in IGV](https://github.com/deLaatLab/peakHiC/raw/master/tutorial/peakHiC_example_igv_snapshot.png)

**Figure 3.** Visualization of peakHiC BigWig tracks containing V4C profiles from 2 example viewpoints

peakHiC loops can also be exported to files / tracks that are compatible with either UCSC or the JuiceBox tool, which is a popular tool for the visualization of Hi-C data. Remember that a plain txt file describing all the loops is already present after running the peakHiC pipeline, which can be further studied with R for instance. To export peakHiC the loops, we choose a PATH (folder) and preFix so that the resulting UCSC Interact BED file and the file with loops in HICCUPS format will be written there.

```{r source}
exportLoops(peakHiCObj,outFilePrefix="/home/geert/localdev/github/peakHiC/RESULTS/Rao_4DN_GM12878_peakHiC_example/rds/loops/peakHiC_example")
```

Below is a screenshot of the exported loops in the file __peakHiC_example_HICCUPS_format.txt__. To obtain this, we downloaded the hg38 assembly version of the GM12878 Hi-C data in juicer (.hic) format from the 4DN, loaded it into Juicebox and imported the file with loops that we exported with peakHiC into Juicebox.

![peakHiC BigWig track in IGV](https://github.com/deLaatLab/peakHiC/raw/master/tutorial/peakHiC_example_Juicebox_snapshot.png)

**Figure 4.** Visualization of peakHiC loops in JuiceBox on top of the GM12878 Hi-C map (hg38)

