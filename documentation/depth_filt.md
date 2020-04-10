# depth_filt()
### Last update: hilldiv 1.5.0
Filter samples based on a minimum sequencing depth.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| countable | A count table (matrix/data.frame) indicating the absolute OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs.
| threshold | A number indicating the minimum sequencing depth required to keep the sample. |

## Examples
````R
data(bat.diet.otutable)

depth_filt(bat.diet.otutable,5000)
depth_filt(bat.diet.otutable,threshold=20000)
````
