# hill_div()
### Last update: hilldiv 1.5.0
Compute neutral or phylogenetic Hill numbers from a single sample (vector) or count table (matrix). Hill numbers or numbers equivalents of diversity indices are diversity measures that compute diversity in effective number of OTUs, i.e. the number of equally abundant OTUs that would be needed to give the same value of diversity.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| count | A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs. |
| qvalue | A positive integer or decimal number (>=0), usually between 0 and 3. |
| tree | An ultrametic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples). |
| dist | A dist object indicating the pairwise distances between samples. THIS FUNCTIONALITY IS NOT IMPLEMENTED YET.  |

## Examples
````R

data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)

#One sample
bat.diet.sample <- bat.diet.otutable[,1]
hill_div(bat.diet.sample,0)
hill_div(bat.diet.sample,qvalue=1)

#One sample (phylogenetic)
names(bat.diet.sample) <- rownames(bat.diet.otutable)
hill_div(bat.diet.sample,1,bat.diet.tree)

#Multiple samples
hill_div(bat.diet.otutable,0)

#Incidence-based
bat.diet.otutable.incidence <- to.incidence(bat.diet.otutable,bat.diet.hierarchy)
hill_div(bat.diet.otutable.incidence,qvalue=1)
hill_div(to.incidence(bat.diet.otutable,bat.diet.hierarchy),1)
````
