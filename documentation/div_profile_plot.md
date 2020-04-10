# div_profile_plot()
### Last update: hilldiv 1.5.0
Create diversity profiles of a single or multiple samples displayed independently or aggregated in groups. Diversity profiles show the relation between the order of diversity (q-value) and the respective Hill numbers, thus providing information about the richness and evenness of a sample at a glance.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| profile | A div.profile() object or a vector/matrix containg diversity profile(s), with columns indicating samples/groups and rows indicating orders of diversity (q-values). |
| colour | A vector of RGB colours e.g. c("#34k235","#99cc00"). The number of vector items, must equal the number of samples or groups that are intended to plot. |
| log | Whether to transform Hill numbers to logarithmic scale (TRUE) or not (FALSE). This is useful when there are large differences between q values (e.g. sharp drop from q=0 to q=1), which might complicate visualization. Default: log="FALSE"  |
| legend | Whether to display the legend (TRUE) or not (FALSE) in diversity profiles containing multiple samples/groups. Default TRUE in multi-sample charts.  |

## Examples
````R

data(bat.diet.otutable)
data(bat.diet.hierarchy)

#One sample example
bat.diet.sample <- bat.diet.otutable[,1]
profile.onesample <- div_profile(count=bat.diet.sample,qvalues=seq(from = 0, to = 5, by = (0.1)))
div_profile_plot(profile.onesample)

#Multiple samples
profile.multiplesamples <- div_profile(bat.diet.otutable)
div_profile_plot(profile.multiplesamples)

#Multiple groups (gamma diversity)
profile.multiplegroups <- div_profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="gamma")
div_profile_plot(profile.multiplegroups)
````

<img align=left src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.one.png" width="350" title="One sample">
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.profile.multiple.png" width="350" title="Multiple samples">
