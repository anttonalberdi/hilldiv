# alpha_div()
### Last update: hilldiv 1.5.0
Broadly speaking, α‐diversity refers to the average diversity of subsystems or samples. It must be highlighted though that the alpha diversity is not obtained by averaging the Hill numbers of the subsystems, but computing the Hill numbers from the averaged basic sums of the subsystems (Chao et al., 2012).

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| otutable | A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs. |
| qvalue | A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals. |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match. |
| weight | A vector indicating the relative weight of each sample. The order needs to be identical to the order of the samples in the OTU table. The values need to sum up to 1. If empty, all samples are weighed the same.  |

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.tree)
alpha_div(countable=bat.diet.otutable,qvalue=1)
alpha_div(countable=bat.diet.otutable,qvalue=1,tree=bat.diet.tree)
weight.vector = rep(1/ncol(bat.diet.otutable),ncol(bat.diet.otutable))
alpha_div(bat.diet.otutable,1,bat.diet.tree,weight.vector)
````
