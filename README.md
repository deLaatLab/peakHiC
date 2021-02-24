# peakHiC - a Hi-C peak calling pipeline

peakHiC is a method written in the R programming language to identify 3D chromatin interactions ("peaks") in Hi-C data. It works by first converting Hi-C contacts into virtual 4C profiles from loci of interest called "viewpoints" and then applying a statistical method to call peaks in these profiles, from replicate Hi-C experiments.

If you have any difficulties using the pipeline, please do not hesitate to contact us (delaatbioinf@hubrecht.eu).

## Citation
Valerio Bianchi, Geert Geeven, Nathan Tucker, Catharina R.E. Hilvering, Amelia W. Hall, Carolina Roselli, Matthew C. Hill, James F. Martin, Kenneth B. Margulies, Patrick T. Ellinor, Wouter de Laat. _Detailed Regulatory Interaction Map of the Human Heart Facilitates Gene Discovery for Cardiovascular Disease._ **bioRxiv 705715**; doi: https://doi.org/10.1101/705715

## Prerequisites

- A Unix like shell (e.g. Bash v3.2+)
- Pairix version 0.3.0+ available from https://github.com/4dn-dcic/pairix
- R v3.5+ available from https://www.r-project.org/.
- The following R packages available from CRAN:
  - argparse
  - config
  - zoo
  - isotone
  - data.table
  - doParallel
- The following R packages available from Bioconductor:
  - GenomicRanges
  - BSgenome of interest (e.g. BSgenome.Hsapiens.UCSC.hg38)
- The peakC package available from https://github.com/deWitLab/peakC/.

## Installation

First choose a folder where to install the pipeline. The R scripts will search for configuration files, HiC data and results relative to this _baseFolder_ of install of peakHiC so it is important to remember this path and set peakHiC up correctly to use it. Here we will use **/home/geert/localdev/github/** as folder where to clone the github repository. After navigating to this folder, we can clone the peakHiC repository using the following command:

```
git clone https://github.com/deLaatLab/peakHiC.git
```

Another path required to setup is the location of the pairix binary. peakHiC uses the 4DN-DCIC tool pairix (see **Prerequisites**) to read HiC data in the pairs format. Below we explain how to configure peakHiC to locate this tool. Please make sure this tool is installed before proceding to the next steps. After installing the pairix tool and the necessary R packages from CRAN and Bioconductor, peakHiC installation should take less than one minute. 

Alternatively, one can use a conda setup as described in the initial comments of the full_run.sh file.

## Configure peakHiC to use example data

Example pairix files for a ~2Mb region on chromosome 1, from the GM12878 HiC dataset published by Rao _et al._ (2014) doi:10.1016/j.cell.2014.11.021 is included to test the peakHiC pipeline. pairix files for other Hi-C datasets are available through the 4DN data portal at https://data.4dnucleome.org/. Run the following commands in the R console to setup peakHiC to use the example data:

```{r source}
baseFolder <- "/home/geert/localdev/github/peakHiC/"
pairixBinary <- "/home/geert/localdev/prog/pairix/bin/pairix"
sourceFile <- paste0(baseFolder,"R/peakHiC_functions.R")
source(sourceFile)
initExampleData(baseFolder=baseFolder,pairixBinary)
```

Please be aware that you need to update the paths in the R code above to point to files and folders on your local machine.

To run peakHiC on a different dataset, adapt the config.yml file and supply the viewpoints of interest and then run the commands commented at the end of the config.yml file.

## Run the pipeline

To run the pipeline, one could either inspect and run the full_run.sh or, alternatively, manually execute the steps as described below:

Navigate in a terminal to the peakHiC _RUN_ folder, so in this example **/home/geert/localdev/github/peakHiC/RUN/** and run the following commands:

```
Rscript createPartitionV4CsByChr.R -chr chr1 -peakHiCObj /home/geert/localdev/github/peakHiC/DATA/example_data/hg38_4DN_Rao_GM12878_peakHiC_example_peakHiCObj.rds
Rscript callPartitionPeaksbyChr.R -chr chr1 -peakHiCObj /home/geert/localdev/github/peakHiC/DATA/example_data/hg38_4DN_Rao_GM12878_peakHiC_example_peakHiCObj.rds
Rscript processLoops.R -chr chr1 -peakHiCObj /home/geert/localdev/github/peakHiC/DATA/example_data/hg38_4DN_Rao_GM12878_peakHiC_example_peakHiCObj.rds
```

Execution time for the above three commands is approximately 200, 560 and 50 seconds respectively using a single thread on a Intel Core i7-6700K CPU at 4.00GHz on a Ubuntu Linux 18.04.3 workstation with 32GB memory. We note that running a single thread should not require more than 2 GB of memory. Running these commands should generate R files with V4C profiles of some example viewpoints and a list of processed loops in **/home/geert/localdev/github/peakHiC/RESULTS/Rao_4DN_GM12878_peakHiC_example/rds/loops/**.

For further analysis and visualization of the example data, please open en read the Tutorial.md (https://github.com/deLaatLab/peakHiC/blob/master/Tutorial.md) in this repository. 