# UqN()

The Jaccard-type overlap quantifies the effective proportion of OTUs or lineages in a system that are shared across all subsystems. Hence, this metric quantifies overlap from the perspective of the overall system. Its corresponding dissimilarity (1 - UqN) quantifies the effective proportion of nonshared OTUs or lineages in the overall system. UqN is integrated in the functions beta.dis() and pair.dis().

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| beta | A beta diversity value based on Hill numbers. |
| qvalue | The q value used to compute the beta diversity. It needs to be a positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals. |
| N | An integer indicating sample size, the number of sampling units to be used to compute the similarity measure.  |

## Examples
````R
UqN(beta=1.24,qvalue=1,N=2)
UqN(1.24,1,2)
````
