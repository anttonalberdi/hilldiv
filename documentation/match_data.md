# match_data()
### Last update: hilldiv 1.5.0
Auxiliary function to match OTUs present in count tables and OTU/ASV phylogenetic trees and ensure phylogenetic Hill numbers are correctly computed.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| countable | A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs. |
| tree | An ultrametic tree of class 'phylo'. |
| output | Whether to output a filtered OTU table (matrix) or a filtered OTU tree (phylo).  |

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.tree)
match_data(bat.diet.otutable,bat.diet.tree,output="countable")
match_data(bat.diet.otutable,bat.diet.tree,output="tree")
````
