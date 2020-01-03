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
  - zoo
  - config
  - data.table
- The following R packages available from Bioconductor:
  - GenomicRanges
  - BSgenome of interest (e.g. BSgenome.Hsapiens.UCSC.hg38)
- The peakC package available from https://github.com/deWitLab/peakC/.


