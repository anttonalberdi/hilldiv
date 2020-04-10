# depth_cov()
### Last update: hilldiv 1.5.0
Coverage of the estimated Hill numbers at different orders of diversity. Assessing whether the sequencing depth of each sample is enough to recover the entire diversity of a sample is an important step to ensure reliable and unbiased comparisons across samples. The function depth.cov() relies on diversity estimations based on Hill numbers (Chao & Jost 2015) to calculate the percentage of estimated diversity covered in each sample.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| abund | A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs. |
| qvalue | A positive integer or decimal number (>=0), usually between 0 and 3. |

## Examples
````R
data(bat.diet.otutable)

depth_cov(bat.diet.otutable,0)
depth_cov(bat.diet.otutable,qvalue=1)
````
