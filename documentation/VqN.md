# VqN()

The Sørensen-type turnover-complement is the complement of the Sørensen‐type turnover, which quantifies the normalized OTU turnover rate with respect to the average subsystem (i.e., alpha), thus provides the proportion of a typical subsystem that changes across subsystems. VqN is integrated in the functions beta.dis() and pair.dis().

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| beta | A beta diversity value based on Hill numbers. |
| N | An integer indicating sample size, the number of sampling units to be used to compute the similarity measure.  |

## Examples
````R
VqN(beta=1.24,N=2)
VqN(1.24,2)
````
