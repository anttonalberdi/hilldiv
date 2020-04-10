# merge_samples()
### Last update: hilldiv 1.5.3
Some diversity analyses, such as comparisons between species that contain multiple samples in the OTU or ASV table, require per-species count information. This function enables merging the count data of samples belonging to the same group, as specified in a hierarchy table.,

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| countable | A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs. |
| hierarchy | A two-column matrix indicating the relation between samples (first column) and groups (second column). |
| relative | Whether to output relative values or not. Default=TRUE. |
| incidence | Whether to transform abundance into incidence data when merging. Default=FALSE.  |

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.hierarchy)
merge_samples(countable=bat.diet.otutable,hierarchy=bat.diet.hierarchy)
merge_samples(bat.diet.otutable,bat.diet.hierarchy)
merge_samples(countable=bat.diet.otutable,hierarchy=bat.diet.hierarchy, incidence=TRUE)
````
