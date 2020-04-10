# beta_dis()
### Last update: hilldiv 1.5.0
Compute dissimilarity or similarity values based on beta diversities (neutral or phylogenetic) and sample size.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| beta | A beta diversity value based on Hill numbers or a div.part() object. |
| qvalue | The q value used to compute the beta diversity. It needs to be a positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals. |
| N | An integer indicating sample size, the number of sampling units to be used to compute the similarity measure.  |
| metric | A character object containing "C", "U", "V" or "S". C: Sørensen‐type overlap or complement. U: Jaccard‐type overlap or complement. V: Sørensen‐type turnover or complement. S: Jaccard‐type turnover or complement. See hilldiv wiki for further information.  |
| type | A character object containing either "similarity" or "dissimilarity". If 'similarity' is used, similarity metrics (0: completely different composition - 1: identical composition) are returned. If 'dissimilarity' is used, dissimilarity metrics (0: identical composition - 1:completely different composition) are returned. |

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.tree)

#Manually indicating beta diversity, order of diversity and sample size
beta_dis(beta=4.5,qvalue=1,N=8)
beta_dis(beta=4.5,qvalue=1,N=8,metric="C",type="similarity")

#Using an object created with the function div_part()
divpartobject <- div_part(bat.diet.otutable,qvalue=0,tree=bat.diet.tree)
beta_dis(divpartobject)
beta_dis(divpartobject,metric="S",type="similarity")
````
