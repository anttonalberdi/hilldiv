# gamma_div()
### Last update: hilldiv 1.5.0
Gamma‐diversity refers to the entire diversity of the system.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| countable | A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs. |
| qvalue | A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals. |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match. |
| weight | A vector indicating the relative weight of each sample. The order needs to be identical to the order of the samples in the OTU table. The values need to sum up to 1. If empty, all samples are weighed the same.  |

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
gamma_div(countable=bat.diet.otutable,qvalue=1)
gamma_div(countable=bat.diet.otutable,qvalue=1,tree=bat.diet.tree)
weight.vector = rep(1/ncol(bat.diet.otutable),ncol(bat.diet.otutable))
gamma_div(bat.diet.otutable,1,bat.diet.tree,weight.vector)
````
