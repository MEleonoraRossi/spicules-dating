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
#         thedefault tail probabilities would have the following equivalent
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
delcol_prior <- 5 # There is one per partition, 5 columns: mu1,mu2,mu3,mu4,mu5

# Path to the directory where the subdirectories where each chain ran for each
# dataset are kept. In this case, we will have the path to the directory where
# the analyses when sampling from the prior took place (i.e., this is the path
# to the `CLK` directory that contains the subdirectories from `1` to `5` as
# we ran 5 chains in this example).
## PRIOR ##
path_prior_NoSpi  <- paste( home_dir, "../01_PAML/01_MCMCtree/analyses/prior/NoSpi/",
                            sep = "" )
path_prior_SpiPor <- paste( home_dir, "../01_PAML/01_MCMCtree/analyses/prior/SpiPor/",
                            sep = "" )

#--------------#
# ANALYSE DATA #
#--------------#

#### ANALYSES NoSpi ####

# 1. Obtain summarise object and generate first convergence plot
prior_sum_NoSpi <- sum_MCMC( num_dirs = num_chains, delcol = delcol_prior,
                             data_dir = path_prior_NoSpi,
                             num_divt = num_divt,
                             node_calib = calib_nodes_all$NoSpi,
                             dataset = "NoSpi",
                             perc = perc, def_samples = def_samples,
                             prior = TRUE,
                             # In this case, I want the output dir to be
                             # created within the same data dir, but it
                             # could be changed!
                             out_dat = path_prior_NoSpi,
                             time_unit = 100 # time unit was scaled by 100
                                             # so that time unit was 100Myr
)
half_chains <- num_chains/2
prior_half1 <- apply( X = prior_sum_NoSpi$mean[1:half_chains,],
                      MARGIN = 2, FUN = mean )
prior_half2 <- apply( X = prior_sum_NoSpi$mean[(half_chains+1):num_chains,],
                      MARGIN = 2, FUN = mean )
pdf( paste( outchecks_dir,
            "ESS_and_chains_convergence/Convergence_plot_prior_NoSpi.pdf",
            sep = "" ),
     paper = "a4" )
plot_convergence( name_dir_dat = "Prior, NoSpi",
                  mean_divt1 = prior_half1,
                  mean_divt2 = prior_half2,
                  num_runs = num_chains )
dev.off()

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
sum_quantiles <- check_quantiles( dat = prior_sum_NoSpi, threshold = th,
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
check_quantiles( dat = prior_sum_NoSpi, threshold = th, main_chain = main_ch,
                 num_chains = num_chains, compare_all = FALSE,
                 outdir = path_prior_NoSpi )
#> END CHECK: no problematic chains!

# 3. Compute ESS with RStan
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
##> more samples may have been collected (e.g., mainly when sampling from the
##> posterior) -- more conservative, but only way the function can be used
##> properly to my knowledge.
ESS_NoSpi <- sum_MCMC_ESS( x = prior_sum_NoSpi$arr4stan,
                           samp_per_chain = prior_sum_NoSpi$samp_per_chain )
# Minimum samples per chain
min( prior_sum_NoSpi$samp_per_chain  )    #20001
# Median samples per chain
median( prior_sum_NoSpi$samp_per_chain  ) #20001
# Maximum samples per chain
max( prior_sum_NoSpi$samp_per_chain  )    #20001
# Number of samples used to compute the ESS
min( prior_sum_NoSpi$samp_per_chain  ) * num_chains #120006
# Show tail-ESS and bulk-ESS
ESS_NoSpi$tab
# Tail-ESS Bulk-ESS
# Med.    30865    33260
# Min.     3359     6695
# Max.    53466    45726
## Calculate Rhat, min and max. Good if max Rhat <= 1.05
min( ESS_NoSpi$stats$Rhat ) # 0.9999812
max( ESS_NoSpi$stats$Rhat ) # 1.001506
# Number of samples collected throughout
dim( prior_sum_NoSpi$all_mcmc ) # 120006     69


#### ANALYSES SpiPor ####

# 1. Obtain summarise object and generate first convergence plot
prior_sum_SpiPor <- sum_MCMC( num_dirs = num_chains, delcol = delcol_prior,
                              data_dir = path_prior_SpiPor,
                              num_divt = num_divt,
                              node_calib = calib_nodes_all$SpiPor,
                              dataset = "SpiPor",
                              perc = perc, def_samples = def_samples,
                              prior = TRUE,
                              # In this case, I want the output dir to be
                              # created within the same data dir, but it
                              # could be changed!
                              out_dat = path_prior_SpiPor,
                              time_unit = 100 # time unit was scaled by 100
                                              # so that time unit was 100 Ma
)
half_chains <- num_chains/2
prior_half1 <- apply( X = prior_sum_SpiPor$mean[1:half_chains,],
                      MARGIN = 2, FUN = mean )
prior_half2 <- apply( X = prior_sum_SpiPor$mean[(half_chains+1):num_chains,],
                      MARGIN = 2, FUN = mean )
pdf( paste( outchecks_dir,
            "ESS_and_chains_convergence/Convergence_plot_prior_SpiPor.pdf",
            sep = "" ),
     paper = "a4" )
plot_convergence( name_dir_dat = "Prior, SpiPor",
                  mean_divt1 = prior_half1,
                  mean_divt2 = prior_half2,
                  num_runs = num_chains )
dev.off()

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
sum_quantiles <- check_quantiles( dat = prior_sum_SpiPor, threshold = th,
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
check_quantiles( dat = prior_sum_SpiPor, threshold = th, main_chain = main_ch,
                 num_chains = num_chains, compare_all = FALSE,
                 outdir = path_prior_SpiPor )
#> END CHECK: no problematic chains!

# 3. Compute ESS with RStan
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
##> more samples may have been collected (e.g., mainly when sampling from the
##> posterior) -- more conservative, but only way the function can be used
##> properly to my knowledge.
ESS_SpiPor <- sum_MCMC_ESS( x = prior_sum_SpiPor$arr4stan,
                           samp_per_chain = prior_sum_SpiPor$samp_per_chain )
# Minimum samples per chain
min( prior_sum_SpiPor$samp_per_chain  )    #20001
# Median samples per chain
median( prior_sum_SpiPor$samp_per_chain  ) #20001
# Maximum samples per chain
max( prior_sum_SpiPor$samp_per_chain  )    #20001
# Number of samples used to compute the ESS
min( prior_sum_SpiPor$samp_per_chain  ) * num_chains #120006
# Show tail-ESS and bulk-ESS
ESS_SpiPor$tab
# Tail-ESS Bulk-ESS
# Tail-ESS Bulk-ESS
# Med.    30662    33336
# Min.     3661     6451
# Max.    59697    58044
## Calculate Rhat, min and max. Good if max Rhat <= 1.05
min( ESS_SpiPor$stats$Rhat ) # 0.9999899
max( ESS_SpiPor$stats$Rhat ) # 1.000565
# Number of samples collected throughout
dim( prior_sum_SpiPor$all_mcmc ) # 120006     69
