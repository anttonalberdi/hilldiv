# pair_dis()
### Last update: hilldiv 1.4.2
Computation of pairwise dissimilarities among samples or groups is a common practice in ecology. The function pair_dis() enables generation of pairwise dissimilarity matrices from count tables, making use of functions div_part() and beta_dis(). It is possible to generate pairwise beta diversity values as well as the four (dis)similarity metrics incorporated in *hilldiv*, i.e. UqN, VqN, SqN and VqN. The function enables generating pairwise distance matrices between samples or any parent group as specified in the hierarchy table.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| countable | A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs. |
| qvalue | A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals. |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match.  |
| hierarchy | A matrix indicating the relation between samples (first column) and parent group(s).  |
| level | Hierarchical level at which to compute the pairwise dissimilarities. Only meaningful if a hierarchy table is provided.  |
| metric | A vector containing any combination of "C", "U", "V" or "S". If not provided, all metrics will be computed.  |

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
pair_dis(bat.diet.otutable,qvalue=1)
pair_dis(bat.diet.otutable,qvalue=1,tree=bat.diet.tree,metric="V")
pair_dis(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy,level="2")

````
