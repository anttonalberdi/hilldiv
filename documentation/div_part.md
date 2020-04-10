# div_part()
### Last update: hilldiv 1.5.0
Diversity partitioning based on Hill numbers. the Hill numbers framework enables the diversity of a system to be partitioned following the multiplicative definition (Jost, 2007; Chao, Chiu, & Hsieh, 2012). Alpha diversity is obtained by computing the Hill numbers from the averaged basic sums of the samples, while gamma diversity is obtained by taking the average of OTU relative abundances across samples, and then computing the Hill numbers of the pooled system. The division of gamma diversity by alpha diversity yields the beta diversity, which quantifies how many times richer an entire system is in effective OTUs (gamma diversity) than its constituent samples are on average (alpha diversity). However, the Hill numbers beta diversity can also be considered an actual diversity value, as the same metric also measures the effective number of equally large and completely distinct samples in a system.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| countable | A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs. |
| qvalue | A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals. |
| hierarchy | A matrix indicating the relation between samples (first column) and parent group(s). |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match.data() if the OTU names do not match.  |

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)

#Two level examples (L1=sample (alpha diversity), L2=whole system (gamma diversity))
div_part(bat.diet.otutable,qvalue=1)
div_part(bat.diet.otutable,qvalue=0,tree=bat.diet.tree)

#Three-level example (L1=sample, L2=species, L3=whole system)
div_part(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
````
