# div_test()
### Last update: hilldiv 1.5.0
Diversity comparison test between groups of samples. The function automatically assesses whether the data meets the properties for parametric statistics and performs the appropriate test accordingly: Students' T, ANOVA, Wilcoxon or Kruskal-Wallis. If the posthoc argument is set as TRUE, multiple group comparisons are complemented with post hoc pairwise tests, either Tukey test (parametric) or Dunn test with Benjamini-Hochberg correction (non-parametric).

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| countable | A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs. |
| qvalue |  A positive integer or decimal number (>=0), usually between 0 and 3. |
| hierarchy | A two-column matrix indicating the relation between samples (first column) and groups (second column). |
| tree | A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function [match.data()](match.data.md) if the OTU names do not match.  |
| posthoc | Whether to run post hoc pairwise analyses or not. If TRUE, an ANOVA will be complemented with a Tukey test and a Kruskal-Wallis test will be complemented with a Dunn test. |

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)

div_test(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
div_test(bat.diet.otutable,qvalue=1,hierarchy=bat.diet.hierarchy,tree=bat.diet.tree)
div_test(bat.diet.otutable,2,bat.diet.hierarchy,bat.diet.tree)
div_test(bat.diet.otutable,qvalue=1,hierarchy=bat.diet.hierarchy,posthoc=TRUE)
````
