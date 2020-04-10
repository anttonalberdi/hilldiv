# is.nested()

Checks whether the hierarchical structure specified in a hierarchy table is nested, which is a requirement for multi-level hierarchical partitioning using the function div.part(). i.e. two samples that belong to a common parent group cannot have different grandparent groups. The best example of nested hierarchy is taxonomy: e.g. two species that belong to the same genus cannot belong to different families.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| hierarchy | A matrix indicating the relation between samples (first column) and parent groups. |

## Examples
````R
data(bat.diet.hierarchy)
is.nested(bat.diet.hierarchy)
````
