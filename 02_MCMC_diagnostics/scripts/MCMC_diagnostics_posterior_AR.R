#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------------------------------#
# LOAD PACKAGES, FUNCTIONS, AND SET ENVIRONMENT #
#-----------------------------------------------#
# Load package needed to automatically find the path to the "scripts"
# directory
library( rstudioapi )
scripts_dir   <- gsub( pattern = "scripts..*", replacement = "scripts/",
                       x = getActiveDocumentContext()$path )
setwd( scripts_dir )
# Load the file with all the functions used throughout this script
source( file = "../../src/Functions.R" )
# Run in-house function to set home directory and output directory for ESS
# and convergence tests
# NOTE: It will create a directory called `plots` and another called
# `ESS_and_chains_convergence` inside the `analyses` directory if you have
# not created them yet
home_dir      <- set_homedir()$home
outchecks_dir <- set_homedir()$ESS
# By now, set the working directory to `home_dir`
setwd( home_dir )

#--------------------------------#
# DEFINE USER'S GLOBAL VARIABLES #
#--------------------------------#
# First, we will define the global variables that match the settings in our 
# analysis.

# 1. Number of chains
num_chains <- 6

# 2. Number of divergence times that have been estimated. One trick to find
# this out quickly is to subtract 1 to the number of species. In this case,
# there are 70 taxa, so the number of internal nodes is `n_taxa-1=70-1=69`.
# Another way to verify this is by opening the `mcmc.txt` file and check the
# header. The first element after `Gen` will have the format of `t_nX`, where
# X will be an integer (i.e., 71). Subtract two to this number (i.e., 71-2=69)
# and this will be your number of divergence times that are parameters of the
# MCMC. Please modify the number below so it fits to the dataset you are using. 
num_divt <- 69
# The following variables will be used to generate traceplots
start_divt <- num_divt + 2
stop_divt  <- start_divt + num_divt - 1

# 3. Number of samples that you specified in the `MCMCtree` control file to 
# collect. NOTE: you may have not collect them all, but do not worry!
def_samples <- 20000

# 4. Quantile percentage that you want to set By default, the variable below is 
# set to 0.975 so the 97.5% and 2.5% quantiles. If you want to change this,
# however, just modify the value.
perc <- 0.975

# 5. Load a semicolon-separated file with info about calibrated nodes. Note that
# each column needs to be separated with semicolons and an extra blank line
# after the last row with calibration information needs to be added (i.e., files
# need to have an extra blank line so R does not complain when reading them). 
# If you add a header, please make sure you name the column elements as 
# `Calib;node;Prior`. If not, the R function below will deal with the header. 
# An example of the format you need to follow to summarise the calibration info
# for each node is the following:
#
# ```
# Calib;node;Prior
# ex_n5;5;ST(5.8300,0.0590,0.1120,109.1240)
# ex_n7;7;B(4.1200,4.5200,0.0250,0.0250)
#
# ```
#
# The first column should have the name of the calibration (e.g., Afrotheria, 
# Laurasiatheria, etc.) as it will help you identify which plot belongs to which
# calibration. The second column is the node used in MCMCtree. The third column
# is the calibration used for that node in MCMCtree format.
# 
# [[ NOTES ABOUT ALLOWED CALIBRATION FORMATS]]
#
# Soft-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and a maximum of 0.8 would with  
#         the default tail probabilities would have the following equivalent 
#         formats:
#         >> B(0.6,0.8) | B(0.6,0.8,0.025,0.025)
#  E.g.2: A calibration with a minimum of 0.6 and a maximum of 0.8 would with  
#         the pL=0.001 and pU=0.025 would have the following format. Note that, 
#         whenever you want to modify either pL or pU, you need to write down 
#         the four parameters in the format of "B(min,max,pL,pU)":
#         >> B(0.6,0.8,0.001,0.025)
#
# Lower-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and the default parameters for
#         p = 0.1, c = 1, pL = 0.025:
#         >> L(0.6) | L(0.6,0.1,1,0.025)
#  E.g.2: A calibration with a hard minimum at 0.6, and so pL = 1e-300. Note 
#         that, whenever you want to modify either pL or pU, you need to write  
#         down the four parameters in the format of "L(min,p,c,pL)":
#         >> L(0.6,0.1,1,1e-300)
#
# Upper-bound calibrations: 
#  E.g.1: A calibration with a maximum of 0.8 and the default parameters for 
#         pU = 0.025:
#         >> U(0.8) | U(0.8,0.025)
#  E.g.2: A calibration with a hard maximum at 0.8, and so pU = 1e-300. Note 
#         that, if you want to modify pU, you need to write down the two
#         parameters in the format of "U(max,pU)":
#         >> U(0.8,1e-300)
#
# ST distributions: 
#  The format accepted has four parameters: xi (location, mean root age), 
#  omega (scale), alpha (shape), nu (df). Accepted format: 
#  >> ST(5.8300,0.0590,0.1120,109.1240)
#
# SN distributions: 
#  The format accepted has three parameters: xi (location, mean root age), 
#  omega (scale), alpha (shape). Accepted format: 
#  >> SN(5.8300,0.0590,0.1120)  
#
#
# The next command executes the `read_calib_f` in-house function, which reads
# your input file (semicolon-separated files). The path to the directory where
# you have saved this file is te argument `main_dir` needs. The argument
# `f_names` requires the file name. Argument `dat` requires a character vector
# to label your dataset. If your input file has a header, please keep
# `head_avail = TRUE`. Otherwise, change this to FALSE.
dat_calibs      <- c( "NoSpi", "SpiPor" )
calib_nodes_all <- read_calib_f( main_dir = paste( home_dir, "calib_files/",
                                                   sep = "" ),
                                 f_names = c( "Calibnodes_NoSpi.csv",
                                              "Calibnodes_SpiPor.csv" ),
                                 dat = dat_calibs,
                                 head_avail = TRUE )

# 5. Number of columns in the `mcmc.txt` that are to be deleted as they do not 
# correspond to sample values for divergence times (i.e., the entries are not 
# names following the format `t_nX`). To figure out this number quickly, you 
# can open the `mcmc.txt` file, read the header, and count the number of `mu*`
# and `sigma2*` elements. Do not count the `lnL` value when looking at 
# `mcmc.txt` files generated when sampling from the posterior -- this is 
# automatically accounted for in the in-house R functions that you will 
# subsequently use. 
# E.g., assuming an MCMC ran under a relaxed-clock model with no partitions, 
# we would see `mu` and `sigma2` columns. Therefore, the variable would be set 
# to `delcol = 2`. Please modify the values below according to your dataset for 
# the `mcmc.txt` file generated when sampling from the prior (`delcol_prior`) 
# and when sampling from the posterior (`delcol_posterior`).
delcol_post <- 10 # There are ten columns: mu1-mu5, sigma2_1-sigma_5

# Path to the directory of each data alignment.
# In this case, we will have the path to the directory where
# the analyses when sampling from the prior took place (i.e., this is the path
# to the `GBM` and `ILN` directories that contains the subdirectories from 
# `1` to `32` as we ran 32 chains).
## POSTERIOR ##
path_NoSpi  <- paste( home_dir, "../01_PAML/01_MCMCtree/analyses/posterior/NoSpi/AR/",
                      sep = "" )
path_SpiPor <- paste( home_dir, "../01_PAML/01_MCMCtree/analyses/posterior/SpiPor/AR/",
                      sep = "" )

#--------------#
# ANALYSE DATA #
#--------------#

#### ANALYSES NoSpi AR ####

# 1. Obtain summarise object and generate first convergence plot
sum_NoSpi <- sum_MCMC( num_dirs = num_chains, delcol = delcol_post,
                       data_dir = path_NoSpi,
                       num_divt = num_divt, node_calib = calib_nodes_all$NoSpi,
                       dataset = "NoSpi_AR",
                       perc = perc, def_samples = def_samples,
                       prior = FALSE,
                       # In this case, I want the output dir to be
                       # created within the same data dir, but it
                       # could be changed!
                       out_dat = path_NoSpi,
                       time_unit = 100 # time unit was scaled by 100
                       # so that time unit was 100 Mya
)
half_chains <- trunc( num_chains/2 )
post_half1  <- apply( X = sum_NoSpi$mean[1:half_chains,],
                      MARGIN = 2, FUN = mean )
post_half2  <- apply( X = sum_NoSpi$mean[(half_chains+1):num_chains,],
                      MARGIN = 2, FUN = mean )
pdf( paste( outchecks_dir,
            "ESS_and_chains_convergence/Convergence_plot_post_NoSpi_AR.pdf",
            sep = "" ),
     paper = "a4" )
plot_convergence( name_dir_dat = "Post - NoSpi+AR",
                  mean_divt1 = post_half1,
                  mean_divt2 = post_half2, num_runs = num_chains )
dev.off()
# Plot traces -- only uncomment and run if enough space in disk!
#
# plot_traces( name_dir_dat = "NoSpi",
#              sum_dat = sum_NoSpi,
#              divt = c(start_divt:stop_divt), n_chains = num_chains,
#              out_dir = path_NoSpi )

#> CHECK: Looking for q97.5% and q2.5% issues across the chains. Basically, this 
#> function compares the divergence times estimated for each node across all 
#> chains. The difference between time estimates is then computed. Differences
#> larger than the threshold set by the user will be flagged, and hence the 
#> corresponding chain will be flagged as problematic.
#> 
#> If you are working with deep phylogenies, you may want to be more generous 
#> with the threshold (e.g., >0.4). Otherwise, a threshold lower than 0.4 can 
#> be too stringent and many chains will be flagged as "problematic". Bear in
#> mind that this threshold should not bee too vague either (e.g., you would
#> accept a different of +-threshold in time estimates at the 97.5% and 2.5%
#> quantile) as the larger the threshold, the larger the differences between
#> the estimates collected across chains for the same nodes.
#> For shallow datasets, you may try a threshold `th = 0.2` so you can better
#> refine which chains you keep. Once you run the ESS checks, you can make sure
#> that Rhat is not larger than 1.05 too.
#> If your analyses has returned flagged chains, please check the output file
#> `check_chains.txt` that will be generated before deciding whether a chain is
#> to be kept or deleted from your analysis.
#> 
#> 1. Set threshold and run preliminary comparison anlaysis across MCMC runs
th <- 0.25
sum_quantiles <- check_quantiles( dat = sum_NoSpi, threshold = th,
                                  num_chains = num_chains, compare_all = TRUE )
#> 2. Find out which is the one that should be used as "main chain" against
#> which the rest should be compared (i.e., the one with less differences in
#> q97.5% and q2.5%)
chain_ind <- which( sum_quantiles$sum_chains[,1] %in% min( sum_quantiles$sum_chains[,1] ) )
main_ch   <- sum_quantiles$sum_chains[chain_ind,2]
if( length( main_ch ) > 1 ){
  # Use by default one of the chains that ran first
  # i.e., starting from chain labelled "1"
  main_ch <- min( main_ch )
}
#> 3. Find out which chains should be removed by using this main chain
#>    In this case, the main chain is chain-32
check_quantiles( dat = sum_NoSpi, threshold = th, main_chain = main_ch,
                 num_chains = num_chains, compare_all = FALSE,
                 outdir = path_NoSpi )
#>
#> END CHECK: There are 3 problematic chains: 1, 3, and 4.

# 2. Get filtered summary object and generate convergence plots with filtered
# dataset
filt_chains <- c( 2, 5, 6 )
write.table( x = t(c( 2, 5, 6 )),
             file = paste( path_NoSpi, "chains_kept.txt", sep = "" ),
             quote = FALSE, row.names = FALSE, col.names = FALSE )
filt_sum <- sum_MCMC( num_dirs = filt_chains, delcol = delcol_post,
                      data_dir = path_NoSpi,
                      num_divt = num_divt, node_calib = calib_nodes_all$NoSpi,
                      dataset = "NoSpi_AR_FILT",
                      perc = perc, def_samples = def_samples,
                      prior = FALSE,
                      # In this case, I want the output dir to be
                      # created within the same data dir, but it
                      # could be changed!
                      out_dat = path_NoSpi,
                      time_unit = 100 # time unit was scaled by 100
                      # so that time unit was 100 Mya
)
num_filt_chains  <- length( filt_chains )
half_chains_filt <- trunc( num_filt_chains/2 ) #1
FILT_half1 <- apply( X = filt_sum$mean[1:2,], MARGIN = 2, FUN = mean )
FILT_half2 <- filt_sum$mean[3,]
pdf( paste( outchecks_dir,
            "ESS_and_chains_convergence/Convergence_plot_post_NoSpi_AR_FILT.pdf",
            sep = ""), 
     paper = "a4" )
plot_convergence( name_dir_dat = "Post FILT - NoSpi+AR",
                  mean_divt1 = FILT_half1,
                  mean_divt2 = FILT_half2, num_runs = num_filt_chains )
dev.off()
# Update objects
sum_NoSpi <- filt_sum
rm( filt_sum )

# 2. Compute ESS with RStan
##> Each column is assumed to be an MCMC. Rows are iterations for parameter X
##> Source explaining why it is preferable than the function in coda:
##> https://nature.berkeley.edu/~pdevalpine/MCMC_comparisons/nimble_MCMC_comparisons.html
##> 
##> We will compute the ESS taking into account all the final filtered chains
##> NOTE: The same number of rows are required to compute the tail-ESS and 
##> bulk-ESS. The second argument of the in-house function `sum_MCMC_ESS` used
##> below is a vector with the number of samples collected for each independent
##> chain. Once the chain with less samples collected has been identified, then
##> we can crop the number of samples in each element of the array to fit 
##> that number. This is essentially what the in-house function `sum_MCMC_ESS`
##> does below. In that way, only the minimum number of samples collected
##> across all independent chains are used for all the chains, even though 
##> more samples have been collected -- more conservative, but only way the 
##> function can be used properly to my knowledge.
ESS_post  <- sum_MCMC_ESS( x = sum_NoSpi$arr4stan,
                           samp_per_chain = sum_NoSpi$samp_per_chain )
# Median samples per chain
median( sum_NoSpi$samp_per_chain  ) # 20001
# Minimum samples per chain
min( sum_NoSpi$samp_per_chain  )    # 20001
# Maximum samples per chain
max( sum_NoSpi$samp_per_chain  )    # 20001
# Number of samples used to compute the ESS
min( sum_NoSpi$samp_per_chain  ) * num_filt_chains # 60003
# Show tail-ESS and bulk-ESS
#>> Both bulk-ESS and tail-ESS should be at least 100 (approximately) per Markov 
#>> Chain in order to be reliable and indicate that estimates of respective
#>> posterior quantiles are reliable.
#>> REF: https://mc-stan.org/rstan/reference/Rhat.html
ESS_post$tab
##      Tail-ESS Bulk-ESS
## Med.     1063      667
## Min.      273      133
## Max.    23413    18287
# Calculate Rhat, min and max. Good if max Rhat <= 1.05
min( ESS_post$stats$Rhat ) # 0.9999731
max( ESS_post$stats$Rhat ) # 1.038014
# Number of samples collected throughout
dim( sum_NoSpi$all_mcmc ) # [1] 60003    69

#### ANALYSES SpiPor AR ####

# 1. Obtain summarise object and generate first convergence plot
sum_SpiPor <- sum_MCMC( num_dirs = num_chains, delcol = delcol_post,
                        data_dir = path_SpiPor,
                        num_divt = num_divt, node_calib = calib_nodes_all$SpiPor,
                        dataset = "SpiPor_AR",
                        perc = perc, def_samples = def_samples,
                        prior = FALSE,
                        # In this case, I want the output dir to be
                        # created within the same data dir, but it
                        # could be changed!
                        out_dat = path_SpiPor,
                        time_unit = 100 # time unit was scaled by 100
                        # so that time unit was 100 Mya
)
half_chains <- num_chains/2
post_half1  <- apply( X = sum_SpiPor$mean[1:half_chains,],
                      MARGIN = 2, FUN = mean )
post_half2  <- apply( X = sum_SpiPor$mean[(half_chains+1):num_chains,],
                      MARGIN = 2, FUN = mean )
pdf( paste( outchecks_dir,
            "ESS_and_chains_convergence/Convergence_plot_post_SpiPor_AR.pdf",
            sep = "" ),
     paper = "a4" )
plot_convergence( name_dir_dat = "Post - SpiPor+AR",
                  mean_divt1 = post_half1,
                  mean_divt2 = post_half2, num_runs = num_chains )
dev.off()
# Plot traces -- only uncomment and run if enough space in disk!
#
# plot_traces( name_dir_dat = "SpiPor",
#              sum_dat = sum_SpiPor,
#              divt = c(start_divt:stop_divt), n_chains = num_chains,
#              out_dir = path_SpiPor )

#> CHECK: Looking for q97.5% and q2.5% issues across the chains. Basically, this 
#> function compares the divergence times estimated for each node across all 
#> chains. The difference between time estimates is then computed. Differences
#> larger than the threshold set by the user will be flagged, and hence the 
#> corresponding chain will be flagged as problematic.
#> 
#> If you are working with deep phylogenies, you may want to be more generous 
#> with the threshold (e.g., >0.4). Otherwise, a threshold lower than 0.4 can 
#> be too stringent and many chains will be flagged as "problematic". Bear in
#> mind that this threshold should not bee too vague either (e.g., you would
#> accept a different of +-threshold in time estimates at the 97.5% and 2.5%
#> quantile) as the larger the threshold, the larger the differences between
#> the estimates collected across chains for the same nodes.
#> For shallow datasets, you may try a threshold `th = 0.2` so you can better
#> refine which chains you keep. Once you run the ESS checks, you can make sure
#> that Rhat is not larger than 1.05 too.
#> If your analyses has returned flagged chains, please check the output file
#> `check_chains.txt` that will be generated before deciding whether a chain is
#> to be kept or deleted from your analysis.
#> 
#> 1. Set threshold and run preliminary comparison anlaysis across MCMC runs
th <- 0.25
sum_quantiles <- check_quantiles( dat = sum_SpiPor, threshold = th,
                                  num_chains = num_chains, compare_all = TRUE )
#> 2. Find out which is the one that should be used as "main chain" against
#> which the rest should be compared (i.e., the one with less differences in
#> q97.5% and q2.5%)
chain_ind <- which( sum_quantiles$sum_chains[,1] %in% min( sum_quantiles$sum_chains[,1] ) )
main_ch   <- sum_quantiles$sum_chains[chain_ind,2]
if( length( main_ch ) > 1 ){
  # Use by default one of the chains that ran first
  # i.e., starting from chain labelled "1"
  main_ch <- min( main_ch )
}
#> 3. Find out which chains should be removed by using this main chain
#>    In this case, the main chain is chain-32
check_quantiles( dat = sum_SpiPor, threshold = th, main_chain = main_ch,
                 num_chains = num_chains, compare_all = FALSE,
                 outdir = path_SpiPor )
#>
#> END CHECK: There are 2 problematic chains: 4 and 5.

# 2. Get filtered summary object and generate convergence plots with filtered
# dataset
filt_chains <- c( 1:3, 6 )
write.table( x = t(c( 1:3, 6 )),
             file = paste( path_SpiPor, "chains_kept.txt", sep = "" ),
             quote = FALSE, row.names = FALSE, col.names = FALSE )
filt_sum <- sum_MCMC( num_dirs = filt_chains, delcol = delcol_post,
                      data_dir = path_SpiPor,
                      num_divt = num_divt, node_calib = calib_nodes_all$SpiPor,
                      dataset = "SpiPor_AR_FILT",
                      perc = perc, def_samples = def_samples,
                      prior = FALSE,
                      # In this case, I want the output dir to be
                      # created within the same data dir, but it
                      # could be changed!
                      out_dat = path_SpiPor,
                      time_unit = 100 # time unit was scaled by 100
                      # so that time unit was 100 Mya
)
num_filt_chains  <- length( filt_chains )
half_chains_filt <- trunc( num_filt_chains/2 ) #2
FILT_half1 <- apply( X = filt_sum$mean[1:2,], MARGIN = 2, FUN = mean )
FILT_half2 <- apply( X = filt_sum$mean[3:4,], MARGIN = 2, FUN = mean )
pdf( paste( outchecks_dir,
            "ESS_and_chains_convergence/Convergence_plot_post_SpiPor_AR_FILT.pdf",
            sep = ""), 
     paper = "a4" )
plot_convergence( name_dir_dat = "Post FILT - SpiPor+AR",
                  mean_divt1 = FILT_half1,
                  mean_divt2 = FILT_half2, num_runs = num_filt_chains )
dev.off()
# Update objects
sum_SpiPor <- filt_sum
rm( filt_sum )

# 2. Compute ESS with RStan
##> Each column is assumed to be an MCMC. Rows are iterations for parameter X
##> Source explaining why it is preferable than the function in coda:
##> https://nature.berkeley.edu/~pdevalpine/MCMC_comparisons/nimble_MCMC_comparisons.html
##> 
##> We will compute the ESS taking into account all the final filtered chains
##> NOTE: The same number of rows are required to compute the tail-ESS and 
##> bulk-ESS. The second argument of the in-house function `sum_MCMC_ESS` used
##> below is a vector with the number of samples collected for each independent
##> chain. Once the chain with less samples collected has been identified, then
##> we can crop the number of samples in each element of the array to fit 
##> that number. This is essentially what the in-house function `sum_MCMC_ESS`
##> does below. In that way, only the minimum number of samples collected
##> across all independent chains are used for all the chains, even though 
##> more samples have been collected -- more conservative, but only way the 
##> function can be used properly to my knowledge.
ESS_post <- sum_MCMC_ESS( x = sum_SpiPor$arr4stan,
                          samp_per_chain = sum_SpiPor$samp_per_chain )
# Median samples per chain
median( sum_SpiPor$samp_per_chain  ) # 20001
# Minimum samples per chain
min( sum_SpiPor$samp_per_chain  )    # 20001
# Maximum samples per chain
max( sum_SpiPor$samp_per_chain  )    # 20001
# Number of samples used to compute the ESS
min( sum_SpiPor$samp_per_chain  ) * num_filt_chains # 80004
# Show tail-ESS and bulk-ESS
#>> Both bulk-ESS and tail-ESS should be at least 100 (approximately) per Markov 
#>> Chain in order to be reliable and indicate that estimates of respective
#>> posterior quantiles are reliable.
#>> REF: https://mc-stan.org/rstan/reference/Rhat.html
ESS_post$tab
##      Tail-ESS Bulk-ESS
## Med.     1348      741
## Min.      411      214
## Max.    33096    29179
# Calculate Rhat, min and max. Good if max Rhat <= 1.05
min( ESS_post$stats$Rhat ) # 1.000045
max( ESS_post$stats$Rhat ) # 1.025086
# Number of samples collected throughout
dim( sum_SpiPor$all_mcmc ) # [1] 80004    69
