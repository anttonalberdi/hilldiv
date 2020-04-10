# pair_dis_plot()
### Last update: hilldiv 1.5.0
The function pair_dis_plot() uses any of the dissimilarity matrices yielded by pair_dis() (e.g. 1-UqN) to visualize it either as a NMDS chart, an associated Shepard plot or a qgraph plot.

## Examples
````R
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)

pairdisres <- pair_dis(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
pair_dis_plot(pairdisres$L2_CqN,type="NMDS")
pair_dis_plot(pairdisres$L1_CqN,hierarchy=bat.diet.hierarchy,type="qgraph")
````
