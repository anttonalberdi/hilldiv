# index_div()
### Last update: hilldiv 1.5.0
Compute neutral (richness, Shannon index, Simpson index) or phylogenetic (Faith's PD, Allen's H, Rao's Q) diversity indices related to Hill numbers from a vector object (one sample) or an OTU table (matrix or data.frame object; multiple samples). A ultrametric tree object (phylo) is necessary to compute phylogenetic diversity indices. Note that if using a tree the tip labels and the 'names' (vectors) or 'rownames' (matrices) need to be identical. Note that if the number of OTUs and samples is high, computing phylodiversities might require considerable time. If the vector or the OTU table columns do not sum to 1, the data is TSS-normalised.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| abund | A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs. |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match. |
| index | Diversity index to be computed ("richness", "shannon", "simpson", "faith", "allen", "rao"). Default without tree argument: index="richness". Default with tree argument: index="faith".  |

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)

#One sample
bat.diet.sample <- bat.diet.otutable[,1]
index_div(bat.diet.sample)
index_div(bat.diet.sample,index="shannon")

#Multiple samples
index_div(bat.diet.otutable)
index_div(bat.diet.otutable,tree=bat.diet.tree,index="faith")

#Incidence-based
bat.diet.otutable.incidence <- to.incidence(bat.diet.otutable,bat.diet.hierarchy)
index_div(bat.diet.otutable.incidence)
index_div(bat.diet.otutable.incidence,index="simpson")
index_div(to.incidence(bat.diet.otutable,bat.diet.hierarchy),tree=bat.diet.tree)
````
