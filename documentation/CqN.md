# CqN

The Sørensen-type overlap quantifies the effective average proportion of a sub‐systems OTUs (or lineages in the case of phylodiversities) that is shared across all subsystems. This is thus a metric that quantifies overlap from the subsystems perspective. Its corresponding dissimilarity measure (1 - CqN) quantifies the effective average proportion of nonshared OTUs or lineages in a system. CqN is integrated in the functions beta.dis() and pair.dis().

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| beta | A beta diversity value based on Hill numbers. |
| qvalue | The q value used to compute the beta diversity. It needs to be a positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals. |
| N | An integer indicating sample size, the number of sampling units to be used to compute the similarity measure.  |

## Examples
````R
CqN(beta=1.24,qvalue=1,N=2)
CqN(1.24,1,2)
````
