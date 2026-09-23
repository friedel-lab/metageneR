# metageneR: Utility Tools to Create and Manipulate Genomic Regions

<!-- badges: start -->
<!-- badges: end -->

The package “metageneR” helps you through your metagene analysis workflow such as
- the creation of genomic regions from GTF, GFF or BED files with options to filter specific features.
- genomic binning of your coverage data with some convenient wrapper functions
- generating highly customizable plots and analyzing the data behind it

## Installation

You can install the development version of metageneR from
[GitHub](https://github.com/) with: (If you cannot use any packages for installing from GitHub, see below)

```r
# install.packages("pak")
pak::pak("friedel-lab/metageneR")
```
or 
```r
# install.packages("devtools")
devtools::install_github("friedel-lab/metageneR")
```
If you need an alternative to use this package as a user- or directory-specific library (e.g. for a remote cluster) without devtools' install_github, pak or similar:

For Linux:
1. Clone the repository locally (where you manage the packages on your own).
2. In RStudio, make sure you are at the repository root.
3. Use ```devtools::build()```. This will generate a tar.gz file.
4. Scp the tar.gz to the remote server.
5. Use ```install.packages("PACKAGENAME.tar.gz", repo = NULL, lib="YOUR_LOCAL_LIB")```
6. Whenever you need to attach the package, use ```library("PACKAGENAME", lib.loc="path/to/YOUR_LOCAL_LIB")```.

## Example

``` r
library(metageneR)

# 1. Annotation: gene bodies with 3 kb upstream and 4 kb downstream windows
make_windows("genes.gtf", upstream = 3000, downstream = 4000,
             filter = filter_gtf_by(type == "gene"),
             path_to_output = "annotations/genes_3kbup_4kbdown.bed")

# 2. Binning: bedgraphs in input_dir named <sample>_pos.bedgraph / <sample>_neg.bedgraph
#    gene body in 200 bins, flanks in 50 fixed bins of 60 bp (up) and 80 bp (down)
bin_genome(input_dir = "bedgraphs/", strand = 1,
           annotation = "annotations/genes_3kbup_4kbdown.bed",
           output_dir = "binning/genes_3kbup_4kbdown",
           bins = 200, fixedBinSizeUpstream = "50:60", fixedBinSizeDownstream = "50:80")

# 3. Plotting: one curve per condition, bin-wise Wilcoxon p-value strips for all pairs
plot_metagene_experiment(plot_dir = "plots", plot_prefix = "genes",
                         run_dir = "binning/genes_3kbup_4kbdown",
                         condition_mapper = data.frame(
                           grep_name    = c("ctrl", "treated"),   # regex matched against file names
                           display_name = c("Control", "Treated")),
                         plot_pairs = "all")
```

## Binning engine and License
Binning is performed by the 'genomeBinner.jar' adapted from binGenome module of the 
Watchdog workflow management system (https://github.com/watchdog-wms/watchdog-wms-modules) [1]. 
The jar is shipped with the package and requires Java 11 or newer. 

metageneR is released under the GPL-3 (see LICENSE.md).


[1] Michael Kluge, Marie-Sophie Friedl, Amrei L Menzel, Caroline C Friedel, Watchdog 2.0: New developments for reusability, reproducibility, and workflow execution, GigaScience, Volume 9, Issue 6, June 2020, giaa068, https://doi.org/10.1093/gigascience/giaa068
