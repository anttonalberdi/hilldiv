# SqN()

The Jaccard-type turnover-complement is thecomplement of the Jaccard‐type turnover, which quantifies the normalized OTU turnover rate with respect to the whole system (i.e. gamma). SqN is integrated in the functions beta.dis() and pair.dis().

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| beta | A beta diversity value based on Hill numbers. |
| N | An integer indicating sample size, the number of sampling units to be used to compute the similarity measure.  |

## Examples
````R
SqN(beta=1.24,N=2)
SqN(1.24,2)
````
