setwd("/Users/anttonalberdi/github/hilldiv/")
library(devtools)
document()
create_package("/Users/anttonalberdi/github/hilldiv/")
build_manual(".")

#Check examples
devtools::check()

##########
# Build for CRAN
##########

#Remove attributes
#xattr -c /Users/jpl786/github/hilldiv/.DS_Store
#chmod -N /Users/jpl786/github/hilldiv/.DS_Store

rm -r /Users/anttonalberdi/Downloads/hilldiv
cp -r /Users/anttonalberdi/github/hilldiv /Users/jpl786/Downloads/
ls -lha /Users/anttonalberdi/Downloads/hilldiv
cd /Users/anttonalberdi/Downloads/
rm -rf hilldiv/.DS_Store
rm -rf hilldiv/.github
rm -rf hilldiv/.git
rm hilldiv/.Rapp.history
rm hilldiv/CHANGELOG.md
rm hilldiv/README.md
rm hilldiv/LICENSE
rm -r hilldiv/documentation
rm -r hilldiv/figures
rm hilldiv/validation.r
ls -lha /Users/anttonalberdi/Downloads/hilldiv

#To check if correct
R CMD check hilldiv

#To compile
R CMD build hilldiv

#Font installation
sudo tlmgr install inconsolata
sudo tlmgr install helvetic

#Inconsolata installation (old) - not necessary
cd /Users/jpl786/Downloads/inconsolata.tds
sudo mkdir -p /usr/local/texlive/2019basic/texmf-dist/web2c
sudo cp -Rfp * /usr/local/texlive/2019basic/texmf-dist
sudo echo Map zi4.map >> /usr/local/texlive/2019basic/texmf-dist/web2c/updmap.cfg
sudo mktexlsr
sudo -H updmap-sys

##########
# INSTALL HILLDIV FROM GITHUB
##########

library(devtools)
install_github("anttonalberdi/hilldiv", force=TRUE)
library(hilldiv)

##########
# TEST ALL FUNCTIONS
##########

####
# alpha_div()
####

data(bat.diet.otutable)
data(bat.diet.tree)
alpha_div(countable=bat.diet.otutable,qvalue=1)
alpha_div(countable=bat.diet.otutable,qvalue=1,tree=bat.diet.tree)
weight.vector = rep(1/ncol(bat.diet.otutable),ncol(bat.diet.otutable))
alpha_div(bat.diet.otutable,1,bat.diet.tree,weight.vector)

####
# beta_dis()
####

data(bat.diet.otutable)
data(bat.diet.tree)

#Manually indicating beta diversity, order of diversity and sample size
beta_dis(beta=4.5,qvalue=1,N=8)
beta_dis(beta=4.5,qvalue=1,N=8,metric="C",type="similarity")

#Using an object created with the function div_part()
divpartobject <- div_part(bat.diet.otutable,qvalue=0,tree=bat.diet.tree)
beta_dis(divpartobject)
beta_dis(divpartobject,metric="S",type="similarity")

####
# copy_filt()
####

data(bat.diet.otutable)

#Remove singletons from all samples
copy_filt(bat.diet.otutable,2)

#Remove OTUs represented by less than 0.01\% of the total reads per sample.
copy_filt(bat.diet.otutable,0.0001)

####
# CqN()
####

CqN(beta=1.24,qvalue=1,N=3)
CqN(1.24,1,3)

####
# depth_cov()
####

data(bat.diet.otutable)
depth_cov(bat.diet.otutable,0)
depth_cov(bat.diet.otutable,qvalue=1)

####
# depth_filt()
####

data(bat.diet.otutable)
depth_filt(bat.diet.otutable,5000)
depth_filt(bat.diet.otutable,threshold=20000)

####
# div_part()
####

data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
#Two level examples (L1=sample (alpha diversity), L2=whole system (gamma diversity))
div_part(bat.diet.otutable,qvalue=1)
div_part(bat.diet.otutable,qvalue=0,tree=bat.diet.tree)
#Three-level example (L1=sample, L2=species, L3=whole system)
div_part(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)

####
# div.profile.plot()
####

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

####
# div_profile()
####

data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
#One sample example
bat.diet.sample <- bat.diet.otutable[,1]
div_profile(count=bat.diet.sample,qvalues=seq(from = 0, to = 5, by = (0.1)))
#One sample example (phylogenetic Hill numbers)
names(bat.diet.sample) <- rownames(bat.diet.otutable)
div_profile(count=bat.diet.sample,qvalues=seq(from = 0, to = 5, by = (0.1)),tree=bat.diet.tree)
#Multiple samples
div_profile(bat.diet.otutable)
#Multiple groups (gamma diversity)
div_profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="gamma")
#Multiple groups (alpha diversity)
div_profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="alpha")

####
# div_test_plot()
####

data(bat.diet.otutable)
data(bat.diet.hierarchy)
divtestres <- div_test(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
div_test_plot(divtestres,chart="box")
div_test_plot(divtestres,chart="violin")
divtest.res.ph <- div_test(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy,posthoc=TRUE)
div_test_plot(divtest.res.ph,chart="jitter",posthoc=TRUE,threshold=0.5)

####
# div_test()
####

data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
div_test(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
div_test(bat.diet.otutable,qvalue=1,hierarchy=bat.diet.hierarchy,tree=bat.diet.tree)
div_test(bat.diet.otutable,2,bat.diet.hierarchy,bat.diet.tree)
div_test(bat.diet.otutable,qvalue=1,hierarchy=bat.diet.hierarchy,posthoc=TRUE)

####
# gamma_div()
####

data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
gamma_div(countable=bat.diet.otutable,qvalue=1)
gamma_div(countable=bat.diet.otutable,qvalue=1,tree=bat.diet.tree)
weight.vector = rep(1/ncol(bat.diet.otutable),ncol(bat.diet.otutable))
gamma_div(bat.diet.otutable,1,bat.diet.tree,weight.vector)

####
# hill_div()
####

data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
#One sample
bat.diet.sample <- bat.diet.otutable[,1]
hill_div(bat.diet.sample,0)
hill_div(bat.diet.sample,qvalue=1)
#One sample (phylogenetic)
names(bat.diet.sample) <- rownames(bat.diet.otutable)
hill_div(bat.diet.sample,1,bat.diet.tree)
#Multiple samples
hill_div(bat.diet.otutable,0)
#Incidence-based
bat.diet.otutable.incidence <- to.incidence(bat.diet.otutable,bat.diet.hierarchy)
hill_div(bat.diet.otutable.incidence,qvalue=1)
hill_div(to.incidence(bat.diet.otutable,bat.diet.hierarchy),1)

####
# index_div()
####

data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
#One sample
bat.diet.sample <- bat.diet.otutable[,1]
index_div(bat.diet.sample)
index_div(bat.diet.sample,index="shannon")
#Multiple samples
index_div(bat.diet.otutable)
index_div(bat.diet.otutable,tree=bat.diet.tree,index="faith")
#Incidence-based
bat.diet.otutable.incidence <- to.incidence(bat.diet.otutable,bat.diet.hierarchy)
index_div(bat.diet.otutable.incidence)
index_div(bat.diet.otutable.incidence,index="simpson")
index_div(to.incidence(bat.diet.otutable,bat.diet.hierarchy),tree=bat.diet.tree)

####
# is.nested()
####

data(bat.diet.hierarchy)
is.nested(bat.diet.hierarchy)

####
# match_data()
####

data(bat.diet.otutable)
data(bat.diet.tree)
match_data(bat.diet.otutable,bat.diet.tree,output="countable")
match_data(bat.diet.otutable,bat.diet.tree,output="tree")

####
# pair_dis_plot()
####

data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
pairdisres <- pair_dis(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy,level="2")
pair_dis_plot(pairdisres$L2_CqN,hierarchy=bat.diet.hierarchy,type="NMDS",level=2)
pair_dis_plot(pairdisres$L2_CqN,hierarchy=bat.diet.hierarchy,type="qgraph",level=2)

####
# pair.dis()
####
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
pair_dis(bat.diet.otutable,qvalue=1)
\dontrun{
pair_dis(bat.diet.otutable,qvalue=1,tree=bat.diet.tree,metric="V")
}
pair_dis(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy,level="2")

####
# SqN()
####

SqN(beta=1.24,N=2)
SqN(1.24,2)

####
# to.incidence()
####

data(bat.diet.otutable)
data(bat.diet.hierarchy)
to.incidence(bat.diet.otutable)
to.incidence(bat.diet.otutable,bat.diet.hierarchy)
to.incidence(bat.diet.otutable,bat.diet.hierarchy,relative=TRUE)
to.incidence(otutable=bat.diet.otutable,hierarchy=bat.diet.hierarchy,relative=TRUE)

####
# tree_depth()
####

data(bat.diet.otutable)
data(bat.diet.tree)
tree_depth(tree=bat.diet.tree,abund=bat.diet.otutable)
tree_depth(bat.diet.tree,bat.diet.otutable)

####
# tss()
####
data(bat.diet.otutable)
tss(bat.diet.otutable)
bat.diet.sample <- bat.diet.otutable[,1]
tss(bat.diet.sample)

####
# UqN()
####

UqN(beta=1.24,qvalue=1,N=2)
UqN(1.24,1,2)

####
# VqN()
####

VqN(beta=1.24,N=2)
VqN(1.24,2)
