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

The R code below will read this file and replace the vpsGR object in peakHiCObj, so that these viewpoints can be analyzed. Based on the genomic partition and restriction fragments defined in the peakHiC object, the __creatVPs__ function will assign a partID and fragID to each viewpoint. 

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

We can also export the peakHiC V4C tracks as BigWig tracks, which can be added to a trackHub for visualization in the UCSC browser or we can directly load them into IGV. 
Below is a screenshot where we loaded the exported BigWig tracks for 2 example VPs into the IGV browser.

![peakHiC BigWig track in IGV](https://github.com/deLaatLab/peakHiC/raw/master/tutorial/peakHiC_example_igv_snapshot.png)

**Figure 1.** Visualization of peakHiC BigWig tracks containing V4C profiles from 2 example viewpoints
