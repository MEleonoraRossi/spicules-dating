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

#-------#
# TASKS #
#-------#
# Read output tree by PhyloBayes (unrooted and with bl)
tt_bl <- ape::read.nexus( "../PhyloBayes_out.nexus" )
# Specify outgroup
outg <- c( "Bolinopsis_ashleyi_pep_cdhit_B", "Nematostella_vectensis_cdhit_N",
           "Pleurobrachia_bachei_pep_cdhit", "Lottia_gigantea_pep_cdhit_LGIG",
           "Homo_sapiens_pep_cdhit_HSAP", "Aurelia_aurita_pep_cdhit_AAUR" )
# Root tree with `ape` R package
tt_root <- ape::root( tt_bl, outgroup = outg, resolve.root = TRUE )
# Output rooted tree
if( ! dir.exists( "reproduce" ) ){
  dir.create( "reproduce" )
}
ape::write.tree( phy = tt_root, file = "reproduce/tree_root_bl_R.tree" )
# Remember to open this file with FigTree and root it again -- this is to make
# sure the tree is rooted as, sometimes, the command above does not do this 
# properly!
