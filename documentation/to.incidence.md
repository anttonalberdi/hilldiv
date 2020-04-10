# to.incidence()

This function transforms an abundance-based count (OTU/ASV) table into an incidence-based vector (if no hierarchy table is specified) or table (if a hierarchy table is specified.)

Although the Hill number framework was originally developed to deal with relative abundance data (i.e., relative number of sequences assigned to each OTU/ASV), it has recently also been applied to incidence data. This means that the relative abundances of the types detected in each of the subsystems (=samples) belonging to a certain system are overlooked, and the diversity of the system is calculated by computing the relative number of detections of a given type across the whole system. Although incidence data are less informative than abundance data, it is both easier to collect, more comparable, and has been extensively used under the niche theory framework

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| otutable | A count table (matrix/data.frame) indicating the absolute or relative OTU abundances of multiple samples. Columns must refer to samples and rows to OTUs. |
| hierarchy | A two-column matrix indicating the relation between samples (first column) and groups (second column).  |
| relative | Whether to transform the incidence vector or matrix to relative (0-1) values. Default: relative=FALSE. |

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.hierarchy)
to.incidence(bat.diet.otutable)
to.incidence(bat.diet.otutable,bat.diet.hierarchy)
to.incidence(bat.diet.otutable,bat.diet.hierarchy,relative=TRUE)
to.incidence(otutable=bat.diet.otutable,hierarchy=bat.diet.hierarchy,relative=TRUE)
````
