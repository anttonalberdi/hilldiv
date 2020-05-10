# hilldiv 1.5.4
### Last update: April 2020

:warning: **This version of hilldiv (1.5.4) is different to the version available at CRAN (1.5.1). The CRAN package will be updated when validating all changes.**

**hilldiv** is an R package that provides a set of functions to assist analysis of diversity for diet reconstruction, microbial community profiling or more general ecosystem characterisation analyses based on Hill numbers, using OTU/ASV tables and associated phylogenetic trees as inputs. The package includes functions for (phylo)diversity measurement, (phylo)diversity profile plotting, (phylo)diversity comparison between samples and groups,  (phylo)diversity partitioning and (dis)similarity measurement. All of these grounded in abundance-based and incidence-based Hill numbers.

The statistical framework developed around Hill numbers encompasses many of the most broadly employed diversity (e.g. richness, Shannon index, Simpson index), phylogenetic diversity (e.g. Faith’s PD, Allen’s H, Rao’s quadratic entropy) and dissimilarity (e.g. Sørensen index, Unifrac distances) metrics. This enables the most common analyses of diversity to be performed while grounded in a single statistical framework. For details about the use of Hill numbers in molecularly characterised biological systems, read the following article.

Alberdi A, Gilbert MTP. (2019). A guide to the application of Hill numbers to DNA‐based diversity analyses. *Molecular Ecology Resources*. 19(4): 804-817. [https://doi.org/10.1111/1755-0998.13014](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13014)

**If you want/need to cite hilldiv:** Alberdi A, Gilbert MTP. 2019. hilldiv: an R package for the integral analysis of diversity based on Hill numbers. bioRxiv, 545665. https://www.biorxiv.org/content/10.1101/545665v1

Some recent papers that used hilldiv:

* Maechler et al. (2020) Decision making and best practices for taxonomy-free eDNA metabarcoding in biomonitoring using Hill numbers. *bioRxiv*. https://www.biorxiv.org/content/10.1101/2020.03.31.017723v1

* Alberdi et al. (2020) DNA metabarcoding and spatial modelling link diet diversification with distribution homogeneity in European bats. *Nature Communications*. https://www.nature.com/articles/s41467-020-14961-2

* Siren et al. (2019) Taxonomic and Functional Characterization of the Microbial Community During Spontaneous in vitro Fermentation of Riesling Must. *Frontiers in Microbiology*. https://doi.org/10.3389/fmicb.2019.00697

* Chang et al. (2019) The roles of morphological traits, resource variation and resource partitioning associated with the dietary niche expansion in the fish‐eating bat *Myotis pilosus*. *Molecular Ecology*. https://doi.org/10.1111/mec.15127

# List of functions
| Function | Short explanation | Documentation |
| ------------- | ------------- | ------------- |
| hill_div()  | Neutral or phylogenetic Hill number computation   | [LINK](documentation/hill_div.md) |
| index_div()  | Neutral or phylogenetic diversity index computation   | [LINK](documentation/index_div.md) |
| div_profile()  | (Phylo)Diversity profiles of individual samples or groups of samples | [LINK](documentation/div_profile.md) |
| div_profile_plot()  | Visual representation of (phylodiversity profiles) | [LINK](documentation/div_profile_plot.md) |
| div_test()  | Diversity comparison between two or multiple groups of samples   | [LINK](documentation/div_test.md) |
| div_test_plot()  | Visual representation of diversity tests | [LINK](documentation/div_test_plot.md) |
| depth_cov()  | Assessment of the sequencing depth per sample   | [LINK](documentation/depth_cov.md) |
| div_part()  | Hierarchical diversity partitioning   | [LINK](documentation/div_part.md) |
| alpha_div()  | Alpha diversity computation   | [LINK](documentation/alpha_div.md) |
| gamma_div()  | Gamma diversity computation   | [LINK](documentation/gamma_div.md) |
| beta_dis()  | (Dis)similarity computation based on beta diversities   | [LINK](documentation/beta_dis.md) |
| pair_dis()  | Pairwise (dis)similarity computation based on beta diversities   | [LINK](documentation/pair_dis.md) |
| pair_dis_plot()  | Visual representation of pairwise (dis)similarities   | [LINK](documentation/pair_dis_plot.md) |
| dis_nmds()  | NMDS representation pairwise (dis)similarities   | To be added |
| UqN()  | Jaccard-type overlap computation from beta diversities | [LINK](documentation/UqN.md) |
| CqN()  | Sørensen-type overlap from beta diversities   | [LINK](documentation/CqN.md) |
| SqN()  | Jaccard-type turnover-complement from beta diversities   | [LINK](documentation/SqN.md) |
| VqN()  | Sørensen-type turnover-complement from beta diversities   | [LINK](documentation/VqN.md) |
| to.incidence()  | Convert an abundance-based count (OTU/ASV) table into an incidence-based object | [LINK](documentation/to.incidence.md) |
| to.occurrences()  | Convert an abundance-based vector or matrix to occurrences | To be added |
| merge_samples()  | Combine samples into groups defined by the hierarchy table. | [LINK](documentation/merge_samples.md) |
| match_data()  | Filter count tables and trees to match OTUs/ASVs in both data files   | [LINK](documentation/match_data.md) |
| depth_filt()  | Filter samples according to a minimum sequencing depth threshold   | [LINK](documentation/depth_filt.md) |
| copy_filt()  | Filter OTUs/ASVs according to a minimum copy number threshold   | [LINK](documentation/copy_filt.md) |
| is.nested()  | Check if hierarchy structure is nested   | [LINK](documentation/is.nested.md) |
| tss()  | Total sum scaling per sample   | [LINK](documentation/tss.md) |
| tree_depth()  | Tree depth computation   | [LINK](documentation/tree_depth.md) |

# Installation
To install **hilldiv** in your R environment, from October 2019, you can rely on the built-in install.packages() function:
````R
install.packages("hilldiv")
library(hilldiv)
````
Or alternatively, you need to 1) install devtools, 2) load devtools library, 3) install **hilldiv** using devtools and 4) finally load **hilldiv** library to your environment.

````R
install.packages("devtools")
library(devtools)
install_github("anttonalberdi/hilldiv")
library(hilldiv)
````
If not installed, it will automatically install the following dependencies: ggplot2, ggpubr, RColorBrewer, data.table, ape, vegan, geiger, qgraph and FSA.

If the console returns the following error:
````R
"tar: Failed to set default locale"
````
Type the following in the console and restart R
````R
system('defaults write org.R-project.R force.LANG en_US.UTF-8')
````

# Updating
If you want to update your local hilldiv library from Github, you can use the following script:
````R
remove.packages("hilldiv")
install_github("anttonalberdi/hilldiv")
library(hilldiv)
````

# Workflow
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/hilldiv_workflow.png" title="Hilldiv workflow">

# Documentation
I am creating extensive documentation about Hill numbers and implementation of *hilldiv* in the [Hilldiv WIKI](https://github.com/anttonalberdi/hilldiv/wiki).

**[0 - Introduction to Hill numbers](https://github.com/anttonalberdi/hilldiv/wiki/0.-Hill-numbers)**: an introduction to Hill numbers. If you don't know the basics about Hill numbers, you should first read this page.

**[1 - Data files](https://github.com/anttonalberdi/hilldiv/wiki/1.-Data-files)**: basic information about the data types implemented in the package hilldiv.

**[2 - Data preprocessing](https://github.com/anttonalberdi/hilldiv/wiki/2.-Data-preprocessing)**: count table, phylogenetic tree and metadata processing for the optimal using the functions [copy_filt()](documentation/copy_filt.md), [depth_cov()](documentation/depth_cov.md), [match_data()](documentation/match_data.md) and [to.incidence()](documentation/to.incidence.md)

**[3.1 - Diversity computation of a single system](https://github.com/anttonalberdi/hilldiv/wiki/3.1-Diversity-computation-of-a-single-system)**: diversity computation of a single system using the functions [hill_div()](documentation/hill_div.md) and [index_div()](documentation/index_div.md).

**[3.2 - Diversity computation and comparison of multiple systems](https://github.com/anttonalberdi/hilldiv/wiki/3.2-Diversity-computation-and-comparison-of-multiple-systems)**: diversity computation and comparison of multiple systems or contrasting groups using the functions [hill_div()](documentation/hill_div.md), [div_test()](documentation/div_test.md) and [div_test_plot()](documentation/div_test_plot.md).

**[3.3 - Diversity profiles](https://github.com/anttonalberdi/hilldiv/wiki/3.3-Diversity-profiles)**: generation of diversity profile tables and plots using functions [div_profile()](documentation/div_profile.md) and [div_profile_plot()](documentation/div_profile_plot.md).

**[4 - Diversity partitioning](https://github.com/anttonalberdi/hilldiv/wiki/4.-Diversity-partitioning)**: hierarchical diversity partitioning using the function [div_part()](documentation/div_part.md).

# References
* Alberdi A., Gilbert M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19(4), 804-817.
* Chao, A. & Jost, L. (2015) Estimating diversity and entropy profiles via discovery rates of new species. Methods in Ecology and Evolution, 6, 873-882.
* Chao, A., Chiu, C.-H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93(9), 2037–2051.
* Chao, A., Chiu, C.-H., & Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society of London. Series B, Biological Sciences, 365(1558), 3599–3609.
* Hill, M. O. (1973). Diversity and Evenness: A Unifying Notation and Its Consequences. Ecology, 54(2), 427–432.
* Hsieh, T. C., & Chao, A. (2017). Rarefaction and Extrapolation: Making Fair Comparison of Abundance-Sensitive Phylogenetic Diversity among Multiple Assemblages. Systematic Biology, 66(1), 100–111.
* Hsieh, T. C., Ma, K. H., & Chao, A. (2016). iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers). Methods in Ecology and Evolution / British Ecological Society, 7(12), 1451–1456.
* Jost, L. (2006). Entropy and diversity. Oikos , 113, 363–375.
* Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88(10), 2427–2439.
* Marcon, E., & Hérault, B. (2015). entropart: An R Package to Measure and Partition Diversity. Journal of Statistical Software, Articles, 67(8), 1–26.
* Tuomisto, H. (2010). A diversity of beta diversities: straightening up a concept gone awry. Part 1. Defining beta diversity as a function of alpha and gamma diversity. Ecography, 33(1), 2–22.
