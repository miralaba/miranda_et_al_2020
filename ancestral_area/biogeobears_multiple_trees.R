#######################################################
# Run BioGeoBEARS on multiple trees
#######################################################

#' @BioGeoBEARS_run_object Set up the inference model you want to run on 
#' each tree, like you would for a normal single-tree run.
#' @newick_fns A list of the Newick files (e.g., you should extract some trees
#' from a BEAST NEXUS MCMC output)
#' @model_name The name you would like added to output filenames
#' @geog_fns A list of corresponding geography files (by default, these are just .geog instead of .newick)
#' @resfns A list of results filenames for .Rdata files, either for saving BioGeoBEARS analyses on each tree, or 
#' for loading previously-saved runs (by default, these are just _model_name.Rdata instead of .newick)
#' @run_or_collect If you want to run BioGeoBEARS on each tree (slow), use "run". If you just want to 
#' collect the results over all the trees, use "collect".  For both, pick "both".
#' @start_treenum Default 1. Change if you want to skip some trees
#' @end_treenum Default is length(newick_fns). Change if you want to run a subset of trees.
#' @runslow If FALSE, old, saved .Rdata results are loaded via load(). Default TRUE generates new .Rdata files and 
#' saves them via save()
#' 


library(optimx)    # (either 2012 or 2013 version, as of January 2014)
library(FD)        # for FD::maxent() (make sure this is up-to-date)
library(snow)      # (if you want to use multicore functionality; prob. better than library(parallel))
#library(parallel)
library(BioGeoBEARS)
library(stringr)
library(cladoRcpp)
library(ggplot2)
library(ggExtra)
library(gridExtra)    # for grid.table() function  


BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$max_range_size = 3
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50 
BioGeoBEARS_run_object$speedup = FALSE
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 3
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
# Set up models see biogeobears.R
check_BioGeoBEARS_run(BioGeoBEARS_run_object)



newick_fns=list.files("newickTREE", full.names = T, recursive = T)
newick_fns=grep(".newick$", newick_fns, value = T)
model_name="CymbDIVALIKE_5AMR3geo3" 
geog_fns="cymb_geog3.data" 
resfns=NULL 
run_or_collect="both" 
start_treenum=1 
end_treenum=length(newick_fns) 
runslow=TRUE 
plot_params=FALSE


# Check the input newick_fns to make sure they end in .newick
num_newick_strings = str_count(string=newick_fns, pattern="\\.newick")
num_newick_strings_equals_1_TF = num_newick_strings == 1
if (sum(num_newick_strings_equals_1_TF) != length(num_newick_strings_equals_1_TF))
{
  error_txt = paste("\n\nERROR in run_bears_optim_on_multiple_trees(): All filenames in 'newick_fns' must have one and only one '.newick'.\nViolators in your 'newick_fns':\n\n", sep="")
  cat(error_txt)
  cat(newick_fns[num_newick_strings_equals_1_TF=FALSE], sep="\n")
  
  stop("Stopping on error in run_bears_optim_on_multiple_trees()")
}
  
# File names to store the state probabilities at nodes and corners
if (model_name == "")
{
  suffix_txt = ""
} else {
  suffix_txt = paste("_", model_name, sep="")
}
stateprobs_at_nodes_fn = paste("state_probs_at_nodes_across_all_trees", suffix_txt, ".Rdata", sep="")
stateprobs_at_corners_fn = paste("state_probs_at_corners_across_all_trees", suffix_txt, ".Rdata", sep="")
  
# Get names of:
# - geography filenames (geog_fns)
# - results filenames (resfns)
#newick_fns = output_trfns
if (is.null(geog_fns))
{
  geog_fns = gsub(pattern="newick", replacement="geog", x=newick_fns)
  cat(geog_fns, ",")
}
if (is.null(resfns))
{
  replacement_txt = paste(suffix_txt, ".Rdata", sep="")
  resfns = gsub(pattern=".newick", replacement=replacement_txt, x=newick_fns)
}
  
if (runslow == TRUE)
{
  row_start = 1
  row_end = 0
  for (i in start_treenum:end_treenum)
    #for (i in 1:length(newick_fns))
    #for (i in 1:2)
    {
      # If you just want to run the inferences and save the results
      if ((run_or_collect == "run") || (run_or_collect == "both"))
      {
        txt = paste("\n\nRunning inference under '", model_name, "' for tree #", i, " of ", start_treenum, "-", end_treenum, "...\n\n", sep="")
        cat(txt)
        
        BioGeoBEARS_run_object$geogfn = geog_fns
        BioGeoBEARS_run_object$trfn = newick_fns[i]
        resfn = resfns[i]
        
        # Check the max number of areas:
        #tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=BioGeoBEARS_run_object$geogfn)
        # max(rowSums(dfnums_to_numeric(tipranges@df), na.rm=TRUE))
        
        res = bears_optim_run(BioGeoBEARS_run_object)
        res    
        
        save(res, file=resfn)
      } # end if (run_or_collect == TRUE)
      
      # Initialize the matrices if it's the first iteration
      if (i == start_treenum)
      {
        # If collecting results, create big empty tables to store the state probabilities at nodes and corners
        if ((run_or_collect == "collect") || (run_or_collect == "both"))
        {
          # For summarizing over the tree collection:
          # Set up an empty table to hold the state
          # probabilities of each run
          # numrows = numtrees * numnodes per tree
          # numcols = numstates
          example_trfn = newick_fns[[start_treenum]]
          example_tr = read.tree(example_trfn)
          numrows = length(example_tr$tip.label) + example_tr$Nnode
          numrows
          # Load results to res
          load(resfns[[1]])
          numstates = ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
          
          # We will make these matrices double-size, just in case trees vary in size
          # (e.g., fossil inclusion or not), then cut the matrices down by excluding 
          # all NAs in the last column.
          
          # Create a big empty matrix (+1 to numcols for the OTUnames)
          # for the node states
          state_probs_at_nodes_across_all_trees = matrix(data=NA, nrow=length(resfns)*numrows*2, ncol=numstates+1)
          dim(state_probs_at_nodes_across_all_trees)
          
          # Create another matrix for the states at the corners (each corner is below
          # a node; thus probably nothing below the root)
          numrows = nrow(res$ML_marginal_prob_each_state_at_branch_bottom_below_node)
          numstates = ncol(res$ML_marginal_prob_each_state_at_branch_bottom_below_node)
          state_probs_at_corners_across_all_trees = matrix(data=NA, nrow=length(resfns)*numrows*2, ncol=numstates+1)
          dim(state_probs_at_corners_across_all_trees)
        } # end if ((run_or_collect == "collect") || (run_or_collect == "both"))
      } # end if (i == start_treenum)
      
      
      
      
      # If you want to run the summary after each tree, or just the summaries
      if ((run_or_collect == "collect") || (run_or_collect == "both"))
      {
        # Processing previously done inferences
        txt = paste("\nProcessing previous inferences under '", model_name, "' for tree #", i, " of ", start_treenum, "-", end_treenum, "...", sep="")
        cat(txt)
        
        
        # Do processing if desired
        # The goal is:
        # For each node and corner, get:
        # 1. The number of times that node appears in the sample
        #    (this should be approximately the PP of the bipartition)
        # 2. Whenever the node does appear, get the state probabilities
        #    at that node.
        # 3. Sum over all state probabilities and divide by #1
        # 4. Result is ancestral state probabilities averaged over the tree
        
        # For each tree, make an sorted text list of the OTUs descending from each node
        # Add to the state probabilities for each node
        trfn = newick_fns[i]
        tr = read.tree(trfn)
        trtable = prt(tr, printflag=FALSE, get_tipnames=TRUE)
        head(trtable)
        
        # Get the BioGeoBEARS results for this tree
        # Loads to "res"
        resfn = resfns[i]
        load(resfn)
        
        names(res)
        state_probs_at_nodes = res$ML_marginal_prob_each_state_at_branch_top_AT_node
        state_probs_at_corners = res$ML_marginal_prob_each_state_at_branch_bottom_below_node
        dim(state_probs_at_corners)
        dim(trtable)
        
        # Store the results
        row_end = row_start-1 + nrow(state_probs_at_nodes)
        rownums = row_start:row_end
        # first columns are the state probabilities
        colnums = 1:numstates
        names_col = numstates + 1
        
        # Store the state probabilities at nodes
        state_probs_at_nodes_across_all_trees[rownums, colnums] = state_probs_at_nodes
        state_probs_at_corners_across_all_trees[rownums, colnums] = state_probs_at_corners
        
        # Store the OTUs descending from each node
        OTUnames = trtable$tipnames
        state_probs_at_nodes_across_all_trees[rownums, names_col] = OTUnames
        state_probs_at_corners_across_all_trees[rownums, names_col] = OTUnames
        
        # update row_start
        row_start = row_end + 1
      } # end if ((run_or_collect == FALSE) or (run_or_collect == "both"))
    } # end for (i in 1:length(newick_fns))
    
# Delete rows that are all NA
allNA <- function(tmprow)
{
  row_is_allNA_TF = all(is.na(tmprow))
  return(row_is_allNA_TF)
}
rows_allNA_TF = apply(X=state_probs_at_nodes_across_all_trees, MARGIN=1, FUN=allNA)
state_probs_at_nodes_across_all_trees = state_probs_at_nodes_across_all_trees[rows_allNA_TF==FALSE,]
state_probs_at_corners_across_all_trees = state_probs_at_corners_across_all_trees[rows_allNA_TF==FALSE,]
    
# Save the states		
txt = paste("\n\nSaving state probabilities at nodes and corners across your tree sample to:\nworking directory: ", getwd(), "\n'", stateprobs_at_nodes_fn, "'\n'", stateprobs_at_corners_fn, "'\n(may be slow, ~1 minute)\n", sep="")
cat(txt)

# After the for-loop, save the ancestral states matrices if you like
save(state_probs_at_nodes_across_all_trees, file=stateprobs_at_nodes_fn)
save(state_probs_at_corners_across_all_trees, file=stateprobs_at_corners_fn)
  

} else {
    
# Load the states
txt = paste("\n\nLoading state probabilities at nodes and corners across your tree sample from:\nworking directory: ", getwd(), "\n'", stateprobs_at_nodes_fn, "'\n'", stateprobs_at_corners_fn, "'\n(may be slow, ~1 minute)\n", sep="")
cat(txt)
# Load to: state_probs_at_nodes_across_all_trees
load(file=stateprobs_at_nodes_fn)
# Load to: state_probs_at_corners_across_all_trees
load(file=stateprobs_at_corners_fn)
  
} # end if runslow==TRUE
  
  
# Also store parameter estimates
optim_results_table = NULL
get_optim_results = TRUE
if (get_optim_results == TRUE)
{
  cat("\n\nGetting ML parameter estimates for trees #", start_treenum, "-", end_treenum, ":\n", sep="")
  for (i in start_treenum:end_treenum)
    #for (i in 1:length(newick_fns))
    #for (i in 1:84)
  {
    cat(i, " ", sep="")
    
    # Get the BioGeoBEARS results for this tree
    # Loads to "res"
    resfn = resfns[i]
    load(resfn)
    
    # Store the parameter estimates
    optim_results_table = rbind(optim_results_table, res$optim_result)
  } # end for (i in 1:length(newick_fns))
} # end if (get_optim_results == TRUE)


optim_results_mean = colMeans(optim_results_table)
optim_results_sd = apply(X=optim_results_table, MARGIN=2, FUN=sd)

# These results are atomic vectors, convert to data.frames as in optimx
tmp_colnames = names(optim_results_mean)
optim_results_mean = data.frame(matrix(data=optim_results_mean, nrow=1))
names(optim_results_mean) = tmp_colnames
optim_results_sd = data.frame(matrix(data=optim_results_sd, nrow=1))
names(optim_results_sd) = tmp_colnames

# Return results
results_on_multiple_trees = list()
  
  
# Add the state probabilities, if you get those...
if ((run_or_collect == "collect") || (run_or_collect == "both"))
{
  results_on_multiple_trees$state_probs_at_nodes_across_all_trees = state_probs_at_nodes_across_all_trees
  results_on_multiple_trees$state_probs_at_corners_across_all_trees = state_probs_at_corners_across_all_trees
}

# Filenames
results_on_multiple_trees$newick_fns = newick_fns
results_on_multiple_trees$geog_fns = geog_fns
results_on_multiple_trees$resfns = resfns

# Optim results (fast)
results_on_multiple_trees$optim_results_table = optim_results_table
results_on_multiple_trees$optim_results_mean = optim_results_mean
results_on_multiple_trees$optim_results_sd = optim_results_sd
save(results_on_multiple_trees, file=paste0("results_on_multiple_trees_", model_name, ".Rdata"))  


# CALCULATE SUMMARY STATISTICS TO COMPARE
# Set up empty tables to hold the statistical results
restable = data.frame(Model=c("DEC", "DECj", "DIVALIKE", "DIVALIKEj", "BAYAREALIKE", "BAYAREALIKEj"),
                      nparam=NA, LnL=NA, sd=NA, d=NA, sd=NA, e=NA, sd=NA, j=NA, sd=NA)

load("results_on_multiple_trees_CymbDEC_G5MR2.Rdata")
resDEC<-results_on_multiple_trees
restable$nparam[1]<-2
restable$LnL[1]<- round(resDEC$optim_results_mean$value, 3)
restable$sd[1] <- round(resDEC$optim_results_sd$value,3)
restable$d[1] <- round(resDEC$optim_results_mean$p1, 3)
restable$sd.1[1] <- round(resDEC$optim_results_sd$p1, 3)
restable$e[1] <- round(resDEC$optim_results_mean$p2, 3)
restable$sd.2[1] <- round(resDEC$optim_results_sd$p2, 3)
restable$j[1] <- 0
restable$sd.3[1] <- 0
DECLnL <- ggplot(resDEC$optim_results_table, aes(value))+geom_histogram(color="black", fill="white")+
          xlab("LnL")+ylab("")+theme_bw(base_family = "serif")
DECparams <- ggplot(resDEC$optim_results_table, aes(p1,p2))+geom_point()+
             xlab("d")+ylab("e")+theme_bw(base_family = "serif")
DECparams <- ggMarginal(DECparams, type = "histogram", color="black", fill="white")



load("results_on_multiple_trees_CymbDECj_G5MR2.Rdata")
resDECj<-results_on_multiple_trees
restable$nparam[2]<-3
restable$LnL[2]<- round(resDECj$optim_results_mean$value, 3)
restable$sd[2] <- round(resDECj$optim_results_sd$value,3)
restable$d[2] <- round(resDECj$optim_results_mean$p1, 3)
restable$sd.1[2] <- round(resDECj$optim_results_sd$p1, 3)
restable$e[2] <- round(resDECj$optim_results_mean$p2, 3)
restable$sd.2[2] <- round(resDECj$optim_results_sd$p2, 3)
restable$j[2] <- round(resDECj$optim_results_mean$p3, 3)
restable$sd.3[2] <- round(resDECj$optim_results_sd$p3, 3)
DECjLnL <- ggplot(resDECj$optim_results_table, aes(value))+geom_histogram(color="black", fill="white")+
          xlab("LnL")+ylab("")+theme_bw(base_family = "serif")
DECjparams <- ggplot(resDECj$optim_results_table, aes(p1,p2))+geom_point()+
             xlab("d")+ylab("e")+theme_bw(base_family = "serif")
DECjparams <- ggMarginal(DECjparams, type = "histogram", color="black", fill="white")




load("results_on_multiple_trees_CymbDIVALIKE_G5MR2.Rdata")
resDIVALIKE<-results_on_multiple_trees
restable$nparam[3]<-2
restable$LnL[3]<- round(resDIVALIKE$optim_results_mean$value, 3)
restable$sd[3] <- round(resDIVALIKE$optim_results_sd$value,3)
restable$d[3] <- round(resDIVALIKE$optim_results_mean$p1, 3)
restable$sd.1[3] <- round(resDIVALIKE$optim_results_sd$p1, 3)
restable$e[3] <- round(resDIVALIKE$optim_results_mean$p2, 3)
restable$sd.2[3] <- round(resDIVALIKE$optim_results_sd$p2, 3)
restable$j[3] <- 0
restable$sd.3[3] <- 0
DIVALIKELnL <- ggplot(resDIVALIKE$optim_results_table, aes(value))+geom_histogram(color="black", fill="white")+
               xlab("LnL")+ylab("")+theme_bw(base_family = "serif")
DIVALIKEparams <- ggplot(resDIVALIKE$optim_results_table, aes(p1,p2))+geom_point()+
                  xlab("d")+ylab("e")+theme_bw(base_family = "serif")
DIVALIKEparams <- ggMarginal(DIVALIKEparams, type = "histogram", color="black", fill="white")



load("results_on_multiple_trees_CymbDIVALIKEj_G5MR2.Rdata")
resDIVALIKEj<-results_on_multiple_trees
restable$nparam[4]<-3
restable$LnL[4]<- round(resDIVALIKEj$optim_results_mean$value, 3)
restable$sd[4] <- round(resDIVALIKEj$optim_results_sd$value,3)
restable$d[4] <- round(resDIVALIKEj$optim_results_mean$p1, 3)
restable$sd.1[4] <- round(resDIVALIKEj$optim_results_sd$p1, 3)
restable$e[4] <- round(resDIVALIKEj$optim_results_mean$p2, 3)
restable$sd.2[4] <- round(resDIVALIKEj$optim_results_sd$p2, 3)
restable$j[4] <- round(resDIVALIKEj$optim_results_mean$p3, 3)
restable$sd.3[4] <- round(resDIVALIKEj$optim_results_sd$p3, 3)
DIVALIKEjLnL <- ggplot(resDIVALIKEj$optim_results_table, aes(value))+geom_histogram(color="black", fill="white")+
                xlab("LnL")+ylab("")+theme_bw(base_family = "serif")
DIVALIKEparams <- ggplot(resDIVALIKEj$optim_results_table, aes(p1,p2))+geom_point()+
                  xlab("d")+ylab("e")+theme_bw(base_family = "serif")
DIVALIKEjparams <- ggMarginal(DIVALIKEjparams, type = "histogram", color="black", fill="white")




load("results_on_multiple_trees_CymbBAYAREALIKE_G5MR2.Rdata")
resBAYAREALIKE<-results_on_multiple_trees
restable$nparam[5]<-2
restable$LnL[5]<- round(resBAYAREALIKE$optim_results_mean$value, 3)
restable$sd[5] <- round(resBAYAREALIKE$optim_results_sd$value,3)
restable$d[5] <- round(resBAYAREALIKE$optim_results_mean$p1, 3)
restable$sd.1[5] <- round(resBAYAREALIKE$optim_results_sd$p1, 3)
restable$e[5] <- round(resBAYAREALIKE$optim_results_mean$p2, 3)
restable$sd.2[5] <- round(resBAYAREALIKE$optim_results_sd$p2, 3)
restable$j[5] <- 0
restable$sd.3[5] <- 0
BAYAREALIKELnL <- ggplot(resBAYAREALIKE$optim_results_table, aes(value))+geom_histogram(color="black", fill="white")+
                  xlab("LnL")+ylab("")+theme_bw(base_family = "serif")
BAYAREALIKEparams <- ggplot(resBAYAREALIKE$optim_results_table, aes(p1,p2))+geom_point()+
                  xlab("d")+ylab("e")+theme_bw(base_family = "serif")
BAYAREALIKEparams <- ggMarginal(BAYAREALIKEparams, type = "histogram", color="black", fill="white")



load("results_on_multiple_trees_CymbBAYAREALIKEj_G5MR2.Rdata")
resBAYAREALIKEj<-results_on_multiple_trees
restable$nparam[6]<-3
restable$LnL[6]<- round(resBAYAREALIKEj$optim_results_mean$value, 3)
restable$sd[6] <- round(resBAYAREALIKEj$optim_results_sd$value,3)
restable$d[6] <- round(resBAYAREALIKEj$optim_results_mean$p1, 3)
restable$sd.1[6] <- round(resBAYAREALIKEj$optim_results_sd$p1, 3)
restable$e[6] <- round(resBAYAREALIKEj$optim_results_mean$p2, 3)
restable$sd.2[6] <- round(resBAYAREALIKEj$optim_results_sd$p2, 3)
restable$j[6] <- round(resBAYAREALIKEj$optim_results_mean$p3, 3)
restable$sd.3[6] <- round(resBAYAREALIKEj$optim_results_sd$p3, 3)
BAYAREALIKEjLnL <- ggplot(resBAYAREALIKEj$optim_results_table, aes(value))+geom_histogram(color="black", fill="white")+
                   xlab("LnL")+ylab("")+theme_bw(base_family = "serif")
BAYAREALIKEparams <- ggplot(resBAYAREALIKEj$optim_results_table, aes(p1,p2))+geom_point()+
                     xlab("d")+ylab("e")+theme_bw(base_family = "serif")
BAYAREALIKEjparams <- ggMarginal(BAYAREALIKEjparams, type = "histogram", color="black", fill="white")



# Model weights of all six models
# With AICcs -- factors in sample size
samplesize <- length(ape::read.tree(resDEC$newick_fns[1])$tip.label)
AICtable <- calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$nparam, samplesize=samplesize)
restable <- cbind(restable, AICtable)
restable <- AkaikeWeights_on_summary_table(restable=restable, colname_to_use<-"AICc")
restable$AICc <- round(restable$AICc, 3)
restable$AICc_wt <- round(restable$AICc_wt, 3)
write.table(restable, "restable_AICc.txt", quote=FALSE, sep="\t")


