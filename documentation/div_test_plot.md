# div_test_plot()
### Last update: hilldiv 1.5.0
Visual comparison between the diversity levels of two or multiple groups of samples, with the option to show the significance values of the pairwise posthoc analyses.

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.hierarchy)

#With no post-hoc analyses
divtestres <- div_test(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
div_test_plot(divtestres,chart="box")
div_test_plot(divtestres,chart="violin")

#With post-hoc analyses
divtest.res.ph <- div_test(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy,posthoc=TRUE)
div_test_plot(divtest.res.ph,chart="jitter",posthoc=TRUE,threshold=0.5)
````
<img align=left src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.test.plot.png" width="400" title="div.test.plot() with pairwise comparisons">
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.test.plot.violin.png" width="400" title="Violin div.test.plot() with pairwise comparisons">
