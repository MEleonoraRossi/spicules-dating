# MCMC diagnostics

## File structure

As soon as our manuscript is accepted for publication, we will update this `README.md` file with a link to download the `analyses` directory with our output files. Nevertheless, you can also run the MCMC diagnostics with your own results if you have decided to run `MCMCtree` on your own with our template files and scripts! Below, you can find the same file structure we used, which you shall require during the MCMC diagnostics if you do not want to update relative/absolute paths that match your own file structure:

```text
01_PAML
 |- 00_CODEML
 |- 01_MCMCtree
     |- analyses
     |    |- posterior
     |    |    |- mcmc_files_[NoSpi|SpiPor]_[AR|IR]  # These directories contain results from MCMC diagnostics
     |    |    |- [NoSpi|SpiPor]
     |    |    |   |- [AR|IR]
     |    |    |       |- [1-6]
     |    |    |       |    |- FigTree.tre
     |    |    |       |    |- in.BV
     |    |    |       |    |- mcmctree*ctl
     |    |    |       |    |- mcmc.txt
     |    |    |       |    |- out.txt
     |    |    |       |    |- SeedUsed
     |    |    |       |- [NoSpi|SpiPor]_[IR|AR]      # These directories contain results from MCMC diagnostics
     |    |    |       |- [NoSpi|SpiPor]_[IR|AR]_FILT # These directories contain results from MCMC diagnostics
     |    |    |       |- *txt                        # These files are output during MCMC diagnostics
     |    |    |- Combine_MCMC.sh # This script is used during MCMC diagnostics
     |    |- prior
     |         |- mcmc_files_[NoSpi|SpiPor] # These directories contain results from MCMC diagnostics
     |         |- [NoSpi|SpiPor]
     |         |   |- [1-6]
     |         |   |     |- FigTree.tre
     |         |   |     |- in.BV
     |         |   |     |- mcmctree*ctl
     |         |   |     |- mcmc.txt
     |         |   |     |- out.txt
     |         |   |     |- SeedUsed
     |         |   |- [NoSpi|SpiPor] # These directories contain results from MCMC diagnostics
     |         |- Combine_MCMC.sh # This script is used during MCMC diagnostics
     |- control_files
     |- in.BV
```

> **NOTE:** you will be able to download the `analyses` directory with the file structure aforementioned once our manuscript has been peer-reviewed.

The file structure that you should have at the beginning of this tutorial inside the `02_MCMC_diagnostics` directory is therefore the following:

```text
02_MCMC_diagnostics
 |- calib_files
 |   |- Calibnodes_NoSpi.csv
 |   |- Calibnodes_SpiPor.csv
 |- dummy_aln
 |   |- dummy_aln.aln
 |- dummy_ctl_files
 |   |- mcmctree_dummy.ctl
 |- plots
 |   |- ESS_and_chains_convergence
 |   |   |- *pdf
 |   |- effVSuser
 |       |- [NoSpi|SpiPor]
 |          |- [NoSpi|SpiPor]/ind/CLK
 |          |    |- CheckUSPvsEP_tn[0-9]*_CLK.png
 |          |- [NoSpi|SpiPor]_[NoSpi|SpiPor]_CLK.jpg
 |- scripts
 |   |- Check_priors_effVSuser.R
 |   |- Combine_MCMC.sh
 |   |- Generate_dummy_aln.R
 |   |- MCMC_diagnostics_posterior_AR.R
 |   |- MCMC_diagnostics_posterior_IR.R
 |   |- MCMC_diagnostics_prior.R
 |- sum_files_post
 |   |- *pdf
 |   |- *tree
 |   |- *tsv
 |- sum_files_prior
 |   |- effVSuser/
 |   |- *pdf
 |   |- *tree
 |   |- *tsv
 |- README.md
```

Once you have the file structures ready, you can follow the instructions below to reproduce our results. Note that you will have output files being generated under both `analyses` and `02_MCMC_diagnostics` directories, and so keeping this file structure is required so that our scripts work! Otherwise, you will have to manually edit them to fix the absolute and relative paths. If you want to reproduce the results we obtained, please run the following commands so that you can generate new output directories and files inside `02_MCMC_diagnostics`:

```sh
# Run from `02_MCMC_diagnostics`
mkdir original_results
mv plots dummy_aln sum_files* original_results/
```

To this end, you shall be able to compare your output with ours if you follow the steps detailed below!

> **NOTE**: we have not given a specific set of commands to place the results we obtained inside directory `analyses` that you are to download as you can decide whether you want to place our results in a specific location prior to running MCMC diagnostics.

## Analysing the samples collected under the prior

We ran an in-house R script, [`MCMC_diagnostics_prior.R`](scripts/MCMC_diagnostics_prior.R) to diagnose the chains that ran in `MCMCtree` when sampling from the prior. In a nutshell, the protocol we followed is the following:

1. Load the `mcmc.txt` files generated after each run.
2. Generate a convergence plot with the unfiltered chains.
3. Find whether there are major differences between the time estimates sampled across the chains for the same nodes in the 97.5% and the 2.5% quantiles. If so, flag and delete said chains.
4. If some chains have not passed the filters mentioned above, create an object with the chains that have passed the filter.
5. Generate a new convergence plot with those chains that passed filters.
6. Calculate the effective sample size for all model parameters and check whether chain convergence has been reached with the chains that have passed filters.

The MCMC diagnostics did not find any of the chains problematic (all details can be found in the R script aforementioned), we could combine the independent chains into a unique samples file to be analysed by the R script [`Check_priors_effVSuser.R`](scripts/Check_priors_effVSuser.R). We used an in-house bash script, [`Combine_MCMC.sh`](scripts/Combine_MCMC.sh), to concatenate in a unique file all the `mcmc.txt` files that correspond to the analysis where chains passed the filters:

```sh
# Run from `02_MCMC_diagnostics/scripts`
cp Combine_MCMC.sh ../../01_PAML/01_MCMCtree/analyses/prior/
cd ../../01_PAML/01_MCMCtree/analyses/prior/
## Variables needed
## arg1   Mame of directory where analyses have taken place (e.g., `fosscal`, `seccal`)
## arg2   Name of the output directory: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3   Filtered chains. E.g., `seq 1 36`, "1 2 5", etc. 
## arg4   Clock model: ILN, GBM, CLK
##
## NOTE: GBM = AR | ILN = IR
./Combine_MCMC.sh NoSpi mcmc_files_NoSpi "`seq 1 6`" CLK
./Combine_MCMC.sh SpiPor mcmc_files_SpiPor "`seq 1 6`" CLK
```

The script above generates two directories called `mcmc_files_NoSpi` and `mcmc_files_SpiPor` inside the `prior` directory, where the `mcmc.txt` with the concatenated samples will be saved for each analysis under the prior. You should see them now too if you are trying to reproduce our analyses.

Then, we generated a dummy alignment with only 2 nucleotides to obtain the summarised timetree with the samples concatenated in the `mcmc.txt` files. In order to do that, we ran our in-house R script [`Generate_dummy_aln.R`](scripts/Generate_dummy_aln.R). If you run it, you can rename the `dummy_aln`you will see a new directory called `dummy_aln` being created, which will contain the dummy alignment we used.

We also generated a dummy control file inside the [`dummy_ctl_files`](dummy_ctl_files) directory with option `print = -1`, which does not start an MCMC but, instead, uses the input files (file with the dummy alignment, calibrated tree file, and concatenated `mcmc.txt` file) to generate a `FigTree.tre` file with the mean estimated divergence times and the corresponding mean CIs using all the samples collected during all the independent MCMCs.

```sh
## NoSpi ##
# Run from `01_MCMCtree/analyses/prior`
base_dir=$( pwd )
cd ../../../../02_MCMC_diagnostics/dummy_ctl_files
ctl_dir=$( pwd )
cd ../../00_data/tree/calib/
tt_dir=$( pwd )
name_tt=`ls *_NoSpi_*tre`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_NoSpi
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
printf "\n"
mv FigTree.tre FigTree_NoSpi.tree
cd $tt_dir
name_tt=`ls *_SpiPori_*tre`
cd $base_dir/mcmc_files_SpiPor
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
printf "\n"
mv FigTree.tre FigTree_SpiPor.tree
cd $base_dir
```

The next step is to plot the user-specified prior (i.e., calibration density) VS the effective prior (i.e., marginal density). We used our in-house R script [`Check_priors_effVSuser.R`](scripts/Check_priors_effVSuser.R) to generate these plots. If you are to run this script with other datasets, however, make sure that your "hard bounds" are not `0.000` in the `Calibnodes_*csv` files and, instead, they are `1e-300` (i.e., while 1e-300 is rounded to `0.000` in the `MCMCtre` output, which can be used to generate the csv files aforementioned, we need `1e-300` to plot distributions in R). To make sure this was not affecting our csv files, we ran the following code snippet:

```sh
# Run from `01_MCMCtree/calib_files`
sed -i 's/0\.0000/1e\-300/g' *csv
```

If you are running the R script [`Check_priors_effVSuser.R`](scripts/Check_priors_effVSuser.R) to reproduce our results, you shall see that a new directory `plots/effVSuser` will have been created. Inside this directory, you will find one directory for each individual dataset with individual plots for each node. In addition, all these plots have been merged into a unique document as well (note: some plots may be too small to see for each node, hence why we also generated individual plots).

Once the MCMC diagnostics finish, it is best to extract the final results when sampling from the prior so that everyone can go through them in a much easier manner:

```sh
# Run from `02_MCMC_diagnostics`
mkdir sum_files_prior
cp -R ../01_PAML/01_MCMCtree/analyses/prior/mcmc_files_*/*tree sum_files_prior/
cp -R ../01_PAML/01_MCMCtree/analyses/prior/*Spi*/*Spi*/*all_mean*tsv sum_files_prior/
cp -R plots/ESS_and_chains_convergence/*prior*pdf sum_files_prior/
cp -R plots/effVSuser sum_files_prior/
```

Now, you can verify that the new output files you have generated match with ours!

## Analysing the samples collected under the posterior

Once we verified that there were no truncation issues or disagreements between the user-specified priors and the corresponding effective priors, we executed `MCMCtree` so that the samples were collected when sampling from the posterior (i.e., using the data) under both the autocorrelated-rates and the independent-rates relaxed-clock models and for both tree hypothesis (i.e., the age of Porifera or Heteroscleromorpha constrained, respectively). You can read more details about the settings in [the `README.md` file inside he `01_MCMCtree` directory](../01_PAML/01_MCMCtree/README.md) and in our manuscript.

Once the analyses finished for all analyses, we ran the `MCMC_diagnostic_posterior_[GBM|ILN].R` files inside the [`scripts` directory](scripts) and followed the same protocol as per the analyses under the prior. Our MCMC diagnostics found some chains that had to be filtered out before we could summarise all the chains together. In that way, you shall find convergence plots which file name ends with "FILT". These plots only include the mean divergence times calculated with only those chains that passed the filters.

Then, we used our in-house bash script, [`Combine_MCMC.sh`](scripts/Combine_MCMC.sh), to concatenate in a unique file all the `mcmc.txt` files that correspond to the analysis carried out with chains that passed our filters:

```sh
# Run from `02_MCMC_diagnostics/scripts`
cp Combine_MCMC.sh ../../01_PAML/01_MCMCtree/analyses/posterior
cd ../../01_PAML/01_MCMCtree/analyses/posterior
## Variables needed
## arg1   Name of directory where analyses have taken place (e.g., conc/fosscal/GBM).
## arg2   Name of the output directory: mcmc_files_conc_fosscal_GBM, mcmc_files_conc_fosscal_ILN, etc.
## arg3   Depends on whether some chains were filtered out or not. E.g., `seq 1 36`, "1 2 5", etc.
## arg4   Clock model: ILN, GBM, CLK.
##
## NOTE 1: GBM = AR | ILN = IR
## NOTE 2: Please only run those commands that correspond to the analyses
## you have run!
./Combine_MCMC.sh NoSpi/AR mcmc_files_NoSpi_AR "2 5 6" GBM
./Combine_MCMC.sh SpiPor/AR mcmc_files_SpiPor_AR "1 2 3 6" GBM
./Combine_MCMC.sh NoSpi/IR mcmc_files_NoSpi_IR "3 4 5" ILN
./Combine_MCMC.sh SpiPor/IR mcmc_files_SpiPor_IR "1 2 4 5 6" ILN
```

Once the scripts above have finished, eight new directories called `mcmc_files_*` will be created inside the `posterior` directory. To infer the final timetrees with the mean time estimates using the samples collected by the filtered chains, we need a control file, the calibrated Newick tree, and the dummy alignment we previously generated in section 2 inside this directory:

```sh
# Run from `01_MCMCtree/analyses/prior`
# Please change directories until
# you are there. 
# Depending on the dataset you are to analyse,
# please change the following variables called 
# `dat_set` and `calib` accoridngly:
dat_set=$( echo SpiPor_IR )  # NoSpi_AR  | NoSpi_IR 
                             # SpiPor_AR | SpiPor_IR 
calib=$( echo SpiPori )      # NoSpi | SpiPori
# Now, run the following commands
# to summarise the samples collected
base_dir=$( pwd )
cd ../../../../02_MCMC_diagnostics/dummy_ctl_files
ctl_dir=$( pwd )
cd ../../00_data/tree/calib/
tt_dir=$( pwd )
name_tt=`ls *_$calib"_"*tre`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_$dat_set
printf "[[ Generating tree file for concatenated \"mcmc.txt\" ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_$dat_set.tree
printf "\n"
cd $base_dir
```

Now, once the MCMC diagnostics have finished, you can extract the results generated when sampling from the posterior so that you can go through them in a much easier way:

```sh
# Run from `02_MCMC_diagnostics`
mkdir sum_files_post
cp -R ../01_PAML/01_MCMCtree/analyses/posterior/mcmc_files_*/*tree sum_files_post/
cp -R ../01_PAML/01_MCMCtree/analyses/posterior/*/*/*/*all_mean*tsv sum_files_post/
cp -R plots/ESS_and_chains_convergence/*post*pdf sum_files_post/
```
