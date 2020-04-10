# tss()

Normalise a vector or count matrix to the range of 0-1. Hill numbers require relative (e.g. OTU1: 0.0012, OTU2: 0.0345) rather than absolute (e.g. OTU1: 1332, OTU2: 5678) abundances. All functions within hilldiv integrate a pre-processing function that convert absolute OTU tables into relative. Nevertheless, OTU tables can be normalised in advance using this function.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| abund | A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs. |

## Examples
````R
data(bat.diet.otutable)
tss(bat.diet.otutable)
bat.diet.sample <- bat.diet.otutable[,1]
tss(bat.diet.sample)
````
