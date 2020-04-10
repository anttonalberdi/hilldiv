# copy_filt()
### Last update: hilldiv 1.5.0
It is highly probable that OTUs with very low representation (number of DNA sequence copies) are artifactual, rather than real biological sequences. Therefore, it is recommendable to remove OTUs represented by a number of sequences below a certain threshold. Although the use of absolute thresholds (e.g. remove all sequences with less than 10 copies) is the most commonly employed strategy, ideally thresholds relative to the sequencing depth should be employed (Alberdi et al. 2018) in order to apply a comparable treatment to samples characterised with a different sequencing depth. The function copy.filt() enables applying both absolute (if threshold value is an integer) and relative (if threshold value is a decimal number) copy number thresholds

## Arguments

| Arguments | Description |
| ------------- | ------------- |
| **abund** |  A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs. |
| **threshold** |  An integer or a decimal number indicating the cut-off threshold. If an integer is provided, an absolute threshold is used (same threshold for all samples). If a decimal number is provided a relative copy number threshold is applied (dependent on the sequencing depth of each sample). |

## Examples
````R
data(bat.diet.otutable)

#Remove singletons from all samples
copy_filt(bat.diet.otutable,2)

#Remove OTUs represented by less than 0.01\% of the total reads per sample.
copy_filt(bat.diet.otutable,0.0001)
````
