

setwd("C:/Users/b0040/Documents/cymbilaimus")


#SIMPLE BIOGEOBEARS SETUP
# Install optimx
install.packages("optimx", dependencies=TRUE, repos="http://cran.rstudio.com")

# Also get snow (for parallel processing)
install.packages("snow")

# Install phylobase
install.packages("phylobase", dependencies=TRUE, repos="http://cran.rstudio.com")

# Install BioGeoBEARS from CRAN 0-cloud:
install.packages("BioGeoBEARS", dependencies=TRUE, repos="http://cran.rstudio.com")


library(gridExtra)    # for grid.table() function  
library(optimx)    # (either 2012 or 2013 version, as of January 2014)
library(FD)        # for FD::maxent() (make sure this is up-to-date)
library(snow)      # (if you want to use multicore functionality; prob. better than library(parallel))
#library(parallel)
library(BioGeoBEARS)
#source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
#calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
#calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)    # slight speedup hopefully

### Look at your phylogeny:
trfn <- "cymb_StarBEAST[ed].txt"
cymb_tree <- read.tree(trfn)
#cymb_tree
plot(cymb_tree)
#cymb_tree <- rotate(cymb_tree, node = xxx)
#nodelabels()
title("Cymbilaimus SpeciesTree")
axisPhylo()

### geography file
geogfn = "cymb_geog3.data"
### Look at your geographic range data:
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges
# Maximum range size observed:
#max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size = 3


#######################################################
#######################################################
# DEC AND DEC+J ANALYSIS
#######################################################
# NOTE: The BioGeoBEARS "DEC" model is identical with 
# the Lagrange DEC model, and should return identical
# ML estimates of parameters, and the same 
# log-likelihoods, for the same datasets.
#
# Ancestral state probabilities at nodes will be slightly 
# different, since BioGeoBEARS is reporting the 
# ancestral state probabilities under the global ML
# model, and Lagrange is reporting ancestral state
# probabilities after re-optimizing the likelihood
# after fixing the state at each node. These will 
# be similar, but not identical. See Matzke (2014),
# Systematic Biology, for discussion.
#
# Also see Matzke (2014) for presentation of the 
# DEC+J model.
#######################################################
#######################################################
# Run DEC
#######################################################
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE      # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O'Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change


# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50     # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = FALSE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE        # if FALSE, use optim() instead of optimx()
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)
BioGeoBEARS_run_object$num_cores_to_use = 2
# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC model
# (nothing to do; defaults)


# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "CymbDEC_5AMR3geo3.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

resDEC$outputs@params_table
# modificar valores m?ximos de BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
# se os valores estimados alcancarem o m?ximo
# p.ex.: BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = valor maior que o m?ximo

#######################################################
# Run DEC+J
#######################################################

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)


resfn = "CymbDECj_5AMR3geo3.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}

resDECj$outputs@params_table
# modificar valores m?ximos de BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
# se os valores estimados alcancarem o m?ximo
# p.ex.: BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = valor maior que o m?ximo




#######################################################
# PDF plots
#######################################################
pdffn = "Cymbilaimus_DEC_vs_DECj_5AMR3geo3.pdf"
pdf(pdffn, width=15, height=20)
#par(oma=c(2,0,0,0))
#par(mai=c(0,0,0,0))
#par(mar=c(0,0,0,0))


#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Cymbilaimus"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, 
                                analysis_titletxt, 
                                addl_params=list("j"), 
                                plotwhat="text", 
                                label.offset=1.5, 
                                tipcex=.9, 
                                statecex=1.1, 
                                splitcex=1.1, 
                                titlecex=1, 
                                plotsplits=F, 
                                cornercoords_loc=scriptdir, 
                                include_null_range=T, 
                                tr=NULL, tipranges=NULL)

# Pie chart
plot_BioGeoBEARS_results(results_object, 
                         analysis_titletxt, 
                         addl_params=list("j"), 
                         plotwhat="pie", 
                         label.offset=1.5, 
                         tipcex=.9, 
                         statecex=1.1, 
                         splitcex=1.1, 
                         titlecex=1, 
                         plotsplits=F, 
                         cornercoords_loc=scriptdir, 
                         include_null_range=T, 
                         tr=NULL, tipranges=NULL)

#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on Cymbilaimus"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, 
                                analysis_titletxt, 
                                addl_params=list("j"), 
                                plotwhat="text", 
                                label.offset=1.5, 
                                tipcex=0.9, 
                                statecex=1.1, 
                                splitcex=1.1, 
                                titlecex=1, 
                                plotsplits=F, 
                                cornercoords_loc=scriptdir, 
                                include_null_range=T, 
                                tr=NULL, tipranges=NULL)

# Pie chart
plot_BioGeoBEARS_results(results_object, 
                         analysis_titletxt, 
                         addl_params=list("j"), 
                         plotwhat="pie", 
                         label.offset=1.5, 
                         tipcex=.9, 
                         statecex=1.1, 
                         splitcex=1.1, 
                         titlecex=1, 
                         plotsplits=F, 
                         cornercoords_loc=scriptdir, 
                         include_null_range=T,
                         tr=NULL, tipranges=NULL)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
#######################################################
# DIVALIKE AND DIVALIKE+J ANALYSIS
#######################################################
#######################################################
# NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
# Ronquist (1997)'s parsimony DIVA. It is a likelihood
# interpretation of DIVA, constructed by modelling DIVA's
# processes the way DEC does, but only allowing the 
# processes DIVA allows (widespread vicariance: yes; subset
# sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
#
# DIVALIKE is a likelihood interpretation of parsimony
# DIVA, and it is "like DIVA" -- similar to, but not
# identical to, parsimony DIVA.
#
# I thus now call the model "DIVALIKE", and you should also. ;-)
#######################################################
#######################################################
#######################################################
# Run DIVALIKE
#######################################################
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O'Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change


# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = FALSE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx()
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)
BioGeoBEARS_run_object$num_cores_to_use = 2
# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run


# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "CymbDIVALIKE_5AMR3geo3.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
}

resDIVALIKE$outputs@params_table
# modificar valores m?ximos de BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
# se os valores estimados alcancarem o m?ximo
# p.ex.: BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = valor maior que o m?ximo


#######################################################
# Run DIVALIKE+J
#######################################################

# Set up DIVALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "CymbDIVALIKEj_5AMR3geo3.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKEj = res
}

resDIVALIKEj$outputs@params_table
# modificar valores m?ximos de BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
# se os valores estimados alcancarem o m?ximo
# p.ex.: BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = valor maior que o m?ximo

#######################################################
# PDF plots
#######################################################
pdffn = "Cymbilaimus_DIVALIKE_vs_DIVALIKEj_5AMR23geo3.pdf"
pdf(pdffn, width=15, height=20)


#######################################################
# Plot ancestral states - DIVALIKE
#######################################################
analysis_titletxt ="BioGeoBEARS DIVALIKE on Cymbilaimus"

# Setup
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, 
                                analysis_titletxt, 
                                addl_params=list("j"), 
                                plotwhat="text", 
                                label.offset=1.5, 
                                tipcex=0.9, 
                                statecex=1.1, 
                                splitcex=1.1, 
                                titlecex=1, 
                                plotsplits=F, 
                                cornercoords_loc=scriptdir, 
                                include_null_range=T, 
                                tr=NULL, tipranges=NULL)

# Pie chart
plot_BioGeoBEARS_results(results_object, 
                         analysis_titletxt, 
                         addl_params=list("j"), 
                         plotwhat="pie", 
                         label.offset=1.5, 
                         tipcex=0.9, 
                         statecex=1.1, 
                         splitcex=1.1, 
                         titlecex=1, 
                         plotsplits=F, 
                         cornercoords_loc=scriptdir, 
                         include_null_range=T, 
                         tr=NULL, tipranges=NULL)

#######################################################
# Plot ancestral states - DIVALIKE+J
#######################################################
analysis_titletxt ="BioGeoBEARS DIVALIKE+J on Cymbilaimus"

# Setup
results_object = resDIVALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, 
                                analysis_titletxt, 
                                addl_params=list("j"), 
                                plotwhat="text", 
                                label.offset=1.5, 
                                tipcex=0.9, 
                                statecex=1.1, 
                                splitcex=1.1, 
                                titlecex=1, 
                                plotsplits=F, 
                                cornercoords_loc=scriptdir, 
                                include_null_range=T, 
                                tr=NULL, tipranges=NULL)

# Pie chart
plot_BioGeoBEARS_results(results_object, 
                         analysis_titletxt, 
                         addl_params=list("j"), 
                         plotwhat="pie", 
                         label.offset=1.5, 
                         tipcex=0.9, 
                         statecex=1.1, 
                         splitcex=1.1, 
                         titlecex=1, 
                         plotsplits=F, 
                         cornercoords_loc=scriptdir, 
                         include_null_range=T,
                         tr=NULL, tipranges=NULL)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
#######################################################
# BAYAREALIKE AND BAYAREALIKE+J ANALYSIS
#######################################################
#######################################################
# NOTE: As with DIVA, the BioGeoBEARS BayArea-like model is 
# not identical with the full Bayesian model implemented 
# in the "BayArea" program of Landis et al. (2013). 
#
# Instead, this is a simplified likelihood interpretation
# of the model.  Basically, in BayArea and BioGeoBEARS-BAYAREALIKE, 
# "d" and "e" work like they do in the DEC model of Lagrange 
# (and BioGeoBEARS), and then BayArea's cladogenesis assumption
# (which is that nothing in particular happens at cladogenesis) is 
# replicated by BioGeoBEARS.
#
# This leaves out 3 important things that are in BayArea:
# 1. Distance dependence (you can add this with a distances 
#    matrix + the "x" parameter in BioGeoBEARS, however)
# 2. A correction for disallowing "e" events that drive
#    a species extinct (a null geographic range)
# 3. The neat Bayesian sampling of histories, which allows
#    analyses on large numbers of areas.
#
# The main purpose of having a "BAYAREALIKE" model is 
# to test the importance of the cladogenesis model on 
# particular datasets. Does it help or hurt the data 
# likelihood if there is no special cladogenesis process?
# 
# BAYAREALIKE is a likelihood interpretation of BayArea,
# and it is "like BayArea" -- similar to, but not
# identical to, Bayesian BayArea.
# I thus now call the model "BAYAREALIKE", and you should also. ;-)
#######################################################
#######################################################
#######################################################
# Run BAYAREALIKE
#######################################################
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O'Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change


# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = FALSE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx()
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)
BioGeoBEARS_run_object$num_cores_to_use = 2
# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "CymbBAYAREALIKE_5AMR3geo3.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
}

resBAYAREALIKE$outputs@params_table
# modificar valores m?ximos de BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
# se os valores estimados alcancarem o m?ximo
# p.ex.: BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = valor maior que o m?ximo

#######################################################
# Run BAYAREALIKE+J
#######################################################
# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond 
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
# slightly for each parameter:
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
#
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
#
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "CymbBAYAREALIKEj_5AMR3geo3.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resBAYAREALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKEj = res
}

resBAYAREALIKEj$outputs@params_table
# modificar valores m?ximos de BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
# se os valores estimados alcancarem o m?ximo
# p.ex.: BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = valor maior que o m?ximo

#######################################################
# PDF plots
#######################################################
pdffn = "Cymbilaimus_BAYAREALIKE_vs_BAYAREALIKEj_5AMR3geo3.pdf"
pdf(pdffn, width=15, height=20)

#######################################################
# Plot ancestral states - BAYAREALIKE
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Cymbilaimus"

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, 
                                analysis_titletxt, 
                                addl_params=list("j"), 
                                plotwhat="text", 
                                label.offset=1.5, 
                                tipcex=0.9, 
                                statecex=1.1, 
                                splitcex=1.1, 
                                titlecex=1, 
                                plotsplits=F, 
                                cornercoords_loc=scriptdir, 
                                include_null_range=T, 
                                tr=NULL, tipranges=NULL)

# Pie chart
plot_BioGeoBEARS_results(results_object, 
                         analysis_titletxt, 
                         addl_params=list("j"), 
                         plotwhat="pie", 
                         label.offset=1.5, 
                         tipcex=0.9, 
                         statecex=1.1, 
                         splitcex=1.1, 
                         titlecex=1, 
                         plotsplits=F, 
                         cornercoords_loc=scriptdir, 
                         include_null_range=T,
                         tr=NULL, tipranges=NULL)

#######################################################
# Plot ancestral states - BAYAREALIKE+J
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J on Cymbilaimus"

# Setup
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, 
                                analysis_titletxt, 
                                addl_params=list("j"), 
                                plotwhat="text", 
                                label.offset=1.5, 
                                tipcex=0.9, 
                                statecex=1.1, 
                                splitcex=1.1, 
                                titlecex=1, 
                                plotsplits=F, 
                                cornercoords_loc=scriptdir, 
                                include_null_range=T, 
                                tr=NULL, tipranges=NULL)

# Pie chart
plot_BioGeoBEARS_results(results_object, 
                         analysis_titletxt, 
                         addl_params=list("j"), 
                         plotwhat="pie", 
                         label.offset=1.5, 
                         tipcex=0.9, 
                         statecex=1.1, 
                         splitcex=1.1, 
                         titlecex=1, 
                         plotsplits=F, 
                         cornercoords_loc=scriptdir, 
                         include_null_range=T,
                         tr=NULL, tipranges=NULL)


dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it






#########################################################################
#########################################################################
#########################################################################
#########################################################################
# 
# CALCULATE SUMMARY STATISTICS TO COMPARE
# DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
# 
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#########################################################################
#########################################################################
# REQUIRED READING:
#
# Practical advice / notes / basic principles on statistical model 
#    comparison in general, and in BioGeoBEARS:
# http://phylo.wikidot.com/advice-on-statistical-model-comparison-in-biogeobears
#########################################################################
#########################################################################

# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

#######################################################
# Statistics -- DEC vs. DEC+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams2 = 2
numparams1 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
# confer the same likelihood on the data. See: Brian O'Meara's webpage:
# http://www.brianomeara.info/tutorials/aic
# ...for an intro to LRT, AIC, and AICc

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams2 = 2
numparams1 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DIVALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams2 = 2
numparams1 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#########################################################################
# ASSEMBLE RESULTS TABLES: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
#########################################################################
teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
restable = put_jcol_after_ecol(restable)
restable

# Look at the results!!
restable
teststable

#######################################################
# Save the results tables for later -- check for e.g.
# convergence issues
#######################################################

# Loads to "restable"
save(restable, file="restable_geogfn3_MR2.Rdata")
load(file="restable_geogfn3_MR2.Rdata")

# Loads to "teststable"
save(teststable, file="teststable_geogfn3_MR2.Rdata")
load(file="teststable_geogfn3_MR2.Rdata")

# Also save to text files
write.table(restable, file="restable.txt", quote=FALSE, sep="\t")
write.table(unlist_df(teststable), file="teststable.txt", quote=FALSE, sep="\t")

#######################################################
# Model weights of all six models
#######################################################
restable2 = restable

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike

# With AICcs -- factors in sample size
samplesize = length(cymb_tree$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike

# Also save to text files
write.table(restable_AIC_rellike, file="restable_AIC_rellike.txt", quote=FALSE, sep="\t")
write.table(restable_AICc_rellike, file="restable_AICc_rellike.txt", quote=FALSE, sep="\t")

# Save with nice conditional formatting
write.table(conditional_format_table(restable_AIC_rellike), file="restable_AIC_rellike_formatted.txt", quote=FALSE, sep="\t")
write.table(conditional_format_table(restable_AICc_rellike), file="restable_AICc_rellike_formatted.txt", quote=FALSE, sep="\t")



#######################################################
# Making final Figure from BioGeoBEARS pretty!!!
#######################################################
# Notes:
# 1. Lines staring with "#" are comments
# 2. This is a script kindly made available. 
# 3. In case you have doubts contact rominassbatista@gmail.com
#######################################################
library(cladoRcpp)

cymb_tree <- read.tree("cymb_StarBEAST[ed].txt")
cymb_tree$tip.label <- c("C.sanctaemariae", "Guiana", "Napo", "Imeri", "Base of Andes", "Cental America", "Inambari", "Tapajós", "Rondônia", "Xingu")
# run this for later when you deal with the results
tips=1:length(cymb_tree$tip.label)
nodes=(length(cymb_tree$tip.label)+1):(length(cymb_tree$tip.label)+cymb_tree$Nnode)

# and this is also to deal with the resuls later
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges=order_tipranges_by_tr(tipranges, cymb_tree)
areas=getareas_from_tipranges_object(tipranges)
numstates=numstates_from_numareas(numareas=length(areas)) 
states_list=areas_list_to_states_list_new(areas)


str(resDIVALIKE)
resDIVALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node

write.csv(resDIVALIKEj$ML_marginal_prob_each_state_at_branch_top_AT_node, file = "resDIVALIKEj_ML_marg_prob.txt")
plot(cymb_tree)
axisPhylo()
nodelabels(cex=0.5)
tiplabels(1:length(cymb_tree$tip.label))


#####################
# DEAL WITH RESULTS #
#####################

# which are the most probable states?
#change the object in 211 and run from 209-412
tips#index for tips
nodes#index for nodes


relprobs_matrix<-resDIVALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node
relprobs_matrix
relprobs_matrix_for_internal_states <-relprobs_matrix[nodes,] 

whichmax<-apply(relprobs_matrix[nodes[1]:nodes[length(nodes)],],1,function(x) which(x==max(x)))
whichmax

second_toreplace<-list()
j<-1
for(i in 1:length(whichmax)){ 
  if (length(whichmax[[i]])==2){
    second_toreplace[[j]]<-whichmax[[i]][2]
    j<-j+1
  }
  whichmax[[i]]<-whichmax[[i]][1]
}

whichmax# whichmax contains now only single values, in case there were two states with the same probability you would
        # need to place the repeated values on the seond max value

names_areas<-areas_list_to_states_list_new(areas=areas, maxareas = 3, include_null_range = TRUE, split_ABC = FALSE)
first_max_states<-names_areas[unlist(whichmax)]#Names for the most probable states
unique(first_max_states)

write.csv(first_max_states, "Cymb_resDIVALIKE_first_max_states.csv")


#Second most probable state
secvalue<-relprobs_matrix
secvalue
i=1
j=1


valuemax<-apply(secvalue[nodes[1]:nodes[length(nodes)],],1,function(x) (max(x)))

for(i in (nodes[1]:nodes[length(nodes)])){
  for (j in (1:length(secvalue[1,]))){
    if (secvalue[i,j]==valuemax[i-(nodes[1]-1)]){
      secvalue[i,j]=0
    }
  }
}

#If there were states with the sampe prob. in the most probable values (here are not) 
#then the other state should be placed as the second most probable value, replacing in the secvalue matrix

whichsecondmax<-apply(secvalue[(nodes[1]:nodes[length(nodes)]),],1,function(x) which(x==max(x)))
whichsecondmax
whichsecondmax<-unlist(whichsecondmax)


second_max<-names_areas[whichsecondmax]
second_max
unique(second_max)

write.csv(second_max, "Cymb_resDIVALIKE_second_max_states.csv")


#Third most probable states
Thirdvalue<-secvalue
valuemax<-apply(Thirdvalue[(nodes[1]:nodes[length(nodes)]),],1,function(x) (max(x)))
for(i in ((nodes[1]:nodes[length(nodes)]))){
  for (j in (1:length(Thirdvalue[1,]))){
    if (Thirdvalue[i,j]==valuemax[i-(nodes[1]-1)]){
      Thirdvalue[i,j]=0
    }   
  }
}



whichthirdmax<-apply(Thirdvalue[nodes[1]:nodes[length(nodes)],],1,function(x) which(x==max(x)))
whichthirdmax<-unlist(whichthirdmax)
third_max<-names_areas[whichthirdmax]
write.csv(third_max, "Cymb_resDIVALIKE_third_max.csv")
unique(third_max)
third_max<-unlist(third_max)


#These would be the corresponding states for all first, second and third states
states_123<-c(first_max_states,second_max,third_max)
states_123
unique(states_123)

#Data with the most probable states for each node

first<-list()
j=1
for (i in (nodes[1]-(nodes[1]-1)):(nodes[length(nodes)]-(nodes[1]-1))){
  value<-whichmax[[j]]
  first[[i]]<-relprobs_matrix[nodes[i],value]
  j=j+1
}
first

#Data with the second most probable states for each node

second<-list()
j=1
for (i in (nodes[1]-(nodes[1]-1)):(nodes[length(nodes)]-(nodes[1]-1))){
  valuese<-whichsecondmax[j]
  second[[i]]<-secvalue[nodes[i],valuese]
  j=j+1
}
second

#Data with the third most probable states for each node

Third<-list()

for (i in (nodes[1]-(nodes[1]-1)):(nodes[length(nodes)]-(nodes[1]-1))){
  valueth<-whichthirdmax[[i]] 
  Third[[i]]<-Thirdvalue[nodes[i],valueth]
}

Third<-unlist(Third)

#Data with the rest of the states combined:

fourth<-list()
for (i in (nodes[1]-(nodes[1]-1)):(nodes[length(nodes)]-(nodes[1]-1)))
{
  fourth[[i]]<-(1-(first[[i]]+second[[i]]+Third[[i]]))
}

fourth

#Now make a new matrix with the number of states plus one add. state
#which is a mix from 4th to last probable state)

matrix<-matrix(c(rep(0,(length(Thirdvalue[,1])*(length(Thirdvalue[1,])+1)))), nrow=length(Thirdvalue[,1]), ncol=length(Thirdvalue[1,])+1)
matrix
data_probs<-cbind(unlist(first),unlist(second),unlist(Third),unlist(fourth), unlist(whichmax), unlist(whichsecondmax),unlist(whichthirdmax), rep(length(Thirdvalue[1,])+1, length(nodes)))

#Now we replace values with probability lower than 0.05 by zero, also replace the number of the corresponding states by "2000" (a non existing state)

for (i in (nodes[1]:nodes[length(nodes)])){
  if (data_probs[i-(nodes[1]-1),1]<0.05){
    data_probs[i-(nodes[1]-1),1]=0
    whichmax[i-(nodes[1]-1)]=2000
  }
  if (data_probs[i-(nodes[1]-1),2]<0.05){
    data_probs[i-(nodes[1]-1),2]=0
    whichsecondmax[i-(nodes[1]-1)]=2000
  }
  if (data_probs[i-(nodes[1]-1),3]<0.05){
    data_probs[i-(nodes[1]-1),3]=0
    whichthirdmax[i-(nodes[1]-1)]=2000
  }
}


#Now we make our new vector of colors
colors_matrix_states<-c(names_areas[unlist(whichmax)],names_areas[whichsecondmax],names_areas[whichthirdmax])
colors_matrix_states<-unique(colors_matrix_states)
colors_matrix_states<-unlist(colors_matrix_states[-(which(colors_matrix_states =="NULL"))])#We remove the state NULL corresponding to a fake state number 2000

#Now we edit our data_probs data frame

fourth_b<-list()
for (i in (nodes[1]-(nodes[1]-1)):length(nodes))
{
  fourth_b[[i]]<-(1-(unlist(data_probs[i,1])+unlist(data_probs[i,2])+unlist(data_probs[i,3])))
}

data_probs<-cbind(data_probs[,1],data_probs[,2],data_probs[,3],unlist(fourth_b), unlist(whichmax), unlist(whichsecondmax),unlist(whichthirdmax), rep(length(Thirdvalue[1,])+1, length(nodes)))


#Now we fill out the empty matrix with the corresponding probabilities


for (i in nodes[1]:nodes[length(nodes)]){
  for (j in 1:(length(Thirdvalue[1,])+1)){
    if (data_probs[(i-(nodes[1]-1)),5]==j){
      matrix[i,j]=data_probs[(i-(nodes[1]-1)),1]
    }
    if (data_probs[(i-(nodes[1]-1)),6]==j){
      matrix[i,j]=data_probs[(i-(nodes[1]-1)),2]
    }
    if (data_probs[(i-(nodes[1]-1)),7]==j){
      matrix[i,j]=data_probs[(i-(nodes[1]-1)),3]
    }
    if (data_probs[(i-(nodes[1]-1)),8]==j){
      matrix[i,j]=data_probs[(i-(nodes[1]-1)),4]
    }
  }
}


relprobs_matrix
for (i in 1:length(tips)){
  for (j in 1:(length(Thirdvalue[1,]))){
    matrix[i,j]=relprobs_matrix[i,j]
  }
}

matrix#Now we have a matrix with all states, last one with least prob. satates - first rows correspond
#to tips, following rows correspond to nodes.

# Lets plot #

#############################################################################################################
#############################################################################################################

pdffn = "Cymb_DIVALIKE_final.pdf"
pdf(pdffn, width=15, height=20)
#par(oma=c(2,0,0,0))
#par(mai=c(0,0,0,0))
#par(mar=c(0,0,0,0))


plot(cymb_tree, no.margin=F, edge.width = 2, cex=2.5, adj=0, label.offset=0.095)

areanames=names(tipranges@df)
areanames
include_null_range = TRUE
states_list_0based_index = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=TRUE)
statenames = areas_list_to_states_list_new(areas=areanames, maxareas=max_range_size, include_null_range=TRUE, split_ABC=FALSE)
colors_matrix = get_colors_for_numareas(length(areas))
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index)
colors_list_for_states
colors_list_for_states[[length(Thirdvalue[1,])+1]]="#000000" #We add an extra color (black) for the less probable states
possible_ranges_list_txt = areas_list_to_states_list_new(areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=TRUE)
states_list_0based_index = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=TRUE)

colors_list_for_states[which(possible_ranges_list_txt=="CA")]<-"#7E3A0A" 
colors_list_for_states[which(possible_ranges_list_txt=="NE")]<-"#ac01e6" 
colors_list_for_states[which(possible_ranges_list_txt=="NW")]<-"#f0a609" 
colors_list_for_states[which(possible_ranges_list_txt=="SW")]<-"#fffe01" 
colors_list_for_states[which(possible_ranges_list_txt=="SE")]<-"#0241e8" 
colors_list_for_states[which(possible_ranges_list_txt=="CANE")]<-"#951D78" 
colors_list_for_states[which(possible_ranges_list_txt=="CANW")]<-"#B76F0A"
colors_list_for_states[which(possible_ranges_list_txt=="CASW")]<-"#BF9B06"
colors_list_for_states[which(possible_ranges_list_txt=="CASE")]<-"#403D79"
colors_list_for_states[which(possible_ranges_list_txt=="NENW")]<-"#CE5378"
colors_list_for_states[which(possible_ranges_list_txt=="NESW")]<-"#D57F74"
colors_list_for_states[which(possible_ranges_list_txt=="NESE")]<-"#5721E7"
colors_list_for_states[which(possible_ranges_list_txt=="NWSW")]<-"#F8D105"
colors_list_for_states[which(possible_ranges_list_txt=="NWSE")]<-"#797379"
colors_list_for_states[which(possible_ranges_list_txt=="SWSE")]<-"#819F75"
colors_list_for_states[which(possible_ranges_list_txt=="CANENW")]<-"#B34A53"
colors_list_for_states[which(possible_ranges_list_txt=="CANESW")]<-"#B86850"
colors_list_for_states[which(possible_ranges_list_txt=="CANESE")]<-"#64299D"
colors_list_for_states[which(possible_ranges_list_txt=="CANWSW")]<-"#CF9E07"
colors_list_for_states[which(possible_ranges_list_txt=="CANWSE")]<-"#7B5F54"
colors_list_for_states[which(possible_ranges_list_txt=="CASWSE")]<-"#807D51"
colors_list_for_states[which(possible_ranges_list_txt=="NENWSW")]<-"#DE8C50"
colors_list_for_states[which(possible_ranges_list_txt=="NENWSE")]<-"#8A4D9D"
colors_list_for_states[which(possible_ranges_list_txt=="NESWSE")]<-"#8F6A9A"
colors_list_for_states[which(possible_ranges_list_txt=="NWSWSE")]<-"#A6A151"




MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties="takefirst")
colors_list_for_states_tips<-colors_list_for_states[1:(length(Thirdvalue[1,]))]
cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states_tips, MLstates)
par(font=1)
nodelabels(node=nodes, pie=matrix[nodes[1]:nodes[length(nodes)],], cex=0.9, cex.lab=1, piecol=colors_list_for_states)
#tiplabels(text=MLstates[tips], tip=tips, bg=cols_byNode[tips], cex=0.8, adj=0)
axisPhylo2(side=1, tick=TRUE, las=1, par(cex=1.5), roundlabels=FALSE, minage=0, pos=0.5, padj=0)
#mtext(text="Ma", side=1, line=1, cex=1, par(font=2))

##Make a legend
#
legend(x=0,y=10, 
       legend=unlist(possible_ranges_list_txt[which(possible_ranges_list_txt %in% colors_matrix_states)]), 
       fill=colors_list_for_states[which(possible_ranges_list_txt %in% colors_matrix_states)],
       #title="Cymbilaimus \n BAYAREALIKE +j \n AIC 34.8", 
       xjust=0, cex = 2, #title.adj=0,
       bty="n")


dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it




