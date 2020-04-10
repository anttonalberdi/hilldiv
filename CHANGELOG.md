# Changelog

### v1.5.3 | February 2020
- Bug corrected in pair_dis().
- Bug corrected in div_part(). 
- Added merge_samples() function.
- Added default removal of OTUs/ASVs with no read representation to copy_filt().
- Added Shepard chart to pair_dis_plot().
- Removed the level argument and modified the functioning of the hierarchy table in pair_dis_plot().
- match_data() now outputs the input tree or table if no modification has been done.

### v1.5.2 | November 2019
- Added automatic naming of hierarchy table columns to avoid problems when plotting NMDSs.
- Added tss() function to tree_depth() to avoid errors when inputing absolute abundance tables.

### v1.5.1 | September 2019
- Corrected the bug at pair_dis() that made two samples with identical composition to yield a beta value of NA instead of 1.

### v1.5.0 | September 2019
- Various modifications to adjust hilldiv to CRAN requirements, including multiple function name changes.

### v1.4.2 | September 2019
- Modification of div.part() for multi-level hierarchical diversity partitioning.
- Modification of beta.dis() to accommodate the changes in div.part().
- Modification of pair.dis() to accommodate the changes in div.part() and beta.dis().
- Addition of auxiliary function is.nested().
- Examples of all functions updated.
- div.test.plot() bug corrected.

### v1.4.1 | August 2019
- Implementation of roxygen2.
- Modification of div.profile() and addition of div.profile.plot() function.
- Added the to.incidence() function.
- Added post hoc tests to the div.test() function.

### v1.3 | July 2019
- Added the function index.div() to compute diversity indices related to Hill numbers.
- Added the match.data() function.
- Added the depth.filt() function.
- Added the auxiliary tree.depth() function.
- Added sample weighting option and corrected "object 'weight.L2' not found" error in pair.dis() function.
- Added automatic tss normalisation to alpha.div() and gamma.div() functions.

### v1.2.2 | June 2019
- Added copy.filt() function to filter OTUs according to absolute or relative copy number thresholds.

### v1.2.1 | May 2019
- Option to plot pairwise mean comparison statistical significance values added to pair.dist.plot().
- Added depth.cov() function for assessment of the sequencing depth per sample based on observed and estimated Hill numbers.

### v1.1 | March 2019
- Dependent "geiger" package added to the Description document.
- Magnify issue corrected in pair.dis.plot.r
