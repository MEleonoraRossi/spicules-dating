#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
# Load package to help find the path to this source file 
library(rstudioapi) 
# Get the path of current open file
path_to_file <- getActiveDocumentContext()$path 
# Get working directory path
wd      <- paste( dirname( path_to_file ), "/", sep = "" )
wd.name <- dirname( path_to_file )
# Set wd. Note that the wd will be the same directory where this 
# script is saved.
setwd( wd )


#----------------------------#
# DEFINE GLOBAL VARS BY USER #
#----------------------------#
# Name of output calibrated tree file ready to be used by `MCMCtree`.
# Note that the file name will have the following format
# "<your_selected_out_name>_calib_MCMCtree.tree".
#
# NOTE: Make always sure that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files. This file needs to be
# already in PHYLIP format. Please follow the same format as used in the 
# example tree file provided.
out_names <- c( "Pori_Cauchy_SM_NoSpi_573", "Pori_Cauchy_SM_SpiPori_573")

# Path to your input tree with calibrations. If in the same directory,
# you only need to write the name. If in a different directory, please
# type the absolute or relative path to this file. You need to include the
# flags within square brackets (e.g., [Eumetazoa]) and write them on the node
# that is to be calibrated.
#
# NOTE: Make always sure that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files.
path_trees <- c( "Pori_Cauchy_SM_NoSpi_573_label.tree",
                "Pori_Cauchy_SM_SpiPori_573_label.tree" )

# Path to your input text file that allows you to match the flags you have 
# used to label the nodes that are to be calibrated with the calibration you 
# want to use in `MCMCtree` format. If in the same directory,
# you only need to write the name. If in a different directory, please
# type the absolute or relative path to this file.
# The format of this file meets the following criteria:
#
#   - No header.
#   - One row per calibration.
#   - No spaces at all.
#   - The pipe character, "|", is used to separate the flag you are using on
#     the node you want to calibrate (left) from the `MCMCtree` calibration
#     that you will use to calibrate the node (right).
#   - The `MCMCtree` calibrations are written in `MCMCtree` format and with 
#     single quotation marks. No spaces.
#   - No spaces when writing the name of your flag either.
# 
# E.g.: row in this text file to calibrate node "EUMETAZOA". The flag used in 
# the tree to locate the node Eumetazoa is "EUMETAZOA" and the `MCMCtree`
# calibration is a soft-bound calibration:
#
# ```
# EUMETAZOA|'B(5.610,5.731,0.025,1e-300)'
# ```
#`
# If in doubt, please follow the same format as used in the example text 
# file provided.
path_textconv <- c( "Calib_converter_NoSpi.txt",
                    "Calib_converter_SpiPori.txt" )

#---------------------------------#
# READ TREE AND CALIBRATIONS FILE #
#---------------------------------#
# Read tree and get phylip header
# NOTE: Make always sure that there is at least one blank line at the 
# end of the tree file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files.
phylip_header <- vector( "numeric", length(path_trees) )
tt            <- vector( "list", length(path_trees) )
for( i in 1:length(path_trees) ){
  tt_name    <- path_trees[i]
  tt[[i]]      <- readLines( tt_name )
  phylip_header[i] <- tt[[i]][1]
  tt[[i]]             <- tt[[i]][2]
}

# Read file converter that will be used to rename calibrations
calibrations <- vector( "list", length(path_textconv) )
for( i in 1:length(path_textconv)){
  tmp <- read.table( file = path_textconv[i],
                     stringsAsFactors = F, sep = "|",
                     blank.lines.skip = T )
  colnames( tmp ) <- c( "name", "MCMCtree calib" )
  calibrations[[i]] <- tmp
}


#--------------------------------#
# REPLACE TAGS WITH CALIBRATIONS #
#--------------------------------#
# Replace calibration names with corresponding calibration
for( i in 1:length(path_textconv) ){
  for( j in 1:(dim( calibrations[[i]] )[1]) ){
    tt[[i]] <- gsub( pattern = paste0("\\[",calibrations[[i]][j,1],"\\]"),
                 x = tt[[i]],
                 replacement = paste( "'", calibrations[[i]][j,2], "'", sep = "" ) )
  }
}


#-------------------------------#
# WRITE CALIBRATED TREE IN FILE #
#-------------------------------#
if( ! dir.exists( "reproduce" ) ){
  dir.create( "reproduce" )
}
out_dir <- "reproduce/"
for( i in 1:length(path_textconv) ){
  write( x = phylip_header[i], file = paste( out_dir, out_names[i], ".tre", sep = "" ) )
  write( x = tt[[i]], file = paste( out_dir, out_names[i], ".tre", sep = "" ),
         append = TRUE )
}

