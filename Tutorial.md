# peakHiC - interactive data analysis and export

If you have followed the steps in the README.md Markdown document and have installed the peakHiC pipeline and have run the commands to create V4C profiles and loops for the example data, the resulting files should now be available for further analysis. Here we will briefly describe how to visualize V4C plots from viewpoints of interest in R and how to export loops and V4C BigWig tracks for further analysis using e.g. the UCSC (https://genome.ucsc.edu/) or IGV (https://igv.org/) genome browsers or the Juicebox (https://github.com/aidenlab/Juicebox) tool developed by the lab of Erez Lieberman Aiden. 

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

![peakHiC BigWig track in IGV](https://github.com/deLaatLab/peakHiC/raw/master/tutorial/peakHiC_example_igv_snapshot.png)
