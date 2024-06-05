[![DOI](https://zenodo.org/badge/810862000.svg)](https://zenodo.org/doi/10.5281/zenodo.11488992)

# Independent origins of spicules reconcile the evolutionary history of sponges (Porifera)

## Introduction

In this repository, you will find the data and the various in-house scripts and tools, together with `CODEML` and `MCMCtree` (two programs that are part of the `PAML` software [[Yang 2007](https://pubmed.ncbi.nlm.nih.gov/17483113/)]), that we used to run a **Bayesian analysis for timetree inference** using a **node-dating approach** and an **approximation to the likelihood calculation** implemented in `MCMCtree` ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) to speed up such analysis.

Before we started these analyses, the following data had already been generated:

* [**Molecular alignment**](00_data/aln/conc/Matrix.phy): partitioned alignment with 103,269 sites for 70 taxa. Additional information regarding the data collected and assembled and how `ModelFinder` ([`Kalyaanamoorthy et al., 2018`](https://www.nature.com/articles/nmeth.4285)) was used to infer the partitions of the full molecular alignment, you can always read the corresponding `Methods` section in the main manuscript and the supplementary material. You will find the individual partitions, which were required to calculate the gradient and the Hessian with `CODEML`, inside [`00_data/aln/part`](00_data/aln/part).
* [**Set of calibrations #1**](00_data/tree/calib/add_calibs/Calib_converter_NoSpi.txt): calibrations based on the fossil record in which the age of clade Heteroscleromorpha has been constrained to test the hypothesis of spicule appearance. More information about how they have been established based on the fossil record in the supplementary material of the manuscript. To see how we have incorporated the calibrations in the tree file, please read [the section "Tree calibration" in the `README.md` file you shall find inside the `00_data` directory](00_data/README.md#tree-calibration).
* [**Set of calibrations #2**](00_data/tree/calib/add_calibs/Calib_converter_SpiPori.txt): calibrations based on the fossil record in which the age of clade Porifera has been constrained to test the hypothesis of spicule appearance. More information about how they have been established based on the fossil record in the supplementary material of the manuscript. To see how we have incorporated the calibrations in the tree file, please read [the section "Tree calibration" in the `README.md` file you shall find inside the `00_data` directory](00_data/README.md#tree-calibration).
* [**Phylogeny**](00_data/tree/raw/PhyloBayes_out.nexus): tree topology with branch lengths inferred with `PhyloBayes` ([Lartillot and Philippe, 2004](http://www.atgc-montpellier.fr/download/papers/cat_2004.pdf)). For more information about how this phylogeny was inferred, you can always read the corresponding `Methods` section in the manuscript.

## What can you do with the content of this repository?

You can use the data and instructions provided in this repository to...

* ... understand how we parsed and formatted our input data.
* ... understand how we ran `PAML` software for timetree inference analysis. E.g., specifying substitution models, selecting the most adequate priors according to our dataset, specifying MCMC settings, etc.
* ... understand how we ran MCMC diagnostics to confidently filter the chains that were run under each type of analysis and assessed chain convergence.
* ... obtain the same results we discuss in our article.

## Workflow

The summary of the workflow we have followed to carry out our analyses is the following:

* Parsing and formatting the input data required to run `PAML` software: files with the sequence alignment and a fixed tree topology.
* Inferring the mean evolutionary rate to specify a sensible rate prior.
* Running `PAML` software for timetree inference:
  * Using various in-house pipelines to set up the working environment, the file structure, and the control files required to run `PAML` software.
  * Running `CODEML` to calculate the branch lengths, the gradient, and the Hessian, vectors and matrix required by `MCMCtree` to enable the approximate likelihood calculation.
  * Running `MCMCtree` with the approximate likelihood calculation enabled for timetree inference. We ran the following analyses to assess the impact that (i) different sets of calibrations based on our hypotheses on spicule formation and (ii) relaxed-clock models (i.e., autocorrelated-rates model or geometric Brownian diffusion model [[Thorne et al. 1998](http://www.ncbi.nlm.nih.gov/pubmed/9866200), [Yang and Rannala 2006](http://www.ncbi.nlm.nih.gov/pubmed/16177230)] and independent log-normal rate model [[Rannala and Yang 2007](http://www.ncbi.nlm.nih.gov/pubmed/17558967), [Lemey et al. 2010](http://www.ncbi.nlm.nih.gov/pubmed/20203288)]) can have on species divergence times estimation:
    * **GBM_SpiPor**: autocorrelated-rates model + set of calibrations #1 (age of clade Porifera constrained).
    * **GBM_NoSpi**: autocorrelated-rates model + set of calibrations #1 (age of clade Heteroscleromorpha constrained).
    * **ILN_SpiPor**: independent-rates model + set of calibrations #2 (age of clade Porifera constrained).
    * **ILN_NoSpi**: independent-rates model + set of calibrations #2 (age of clade Heteroscleromorpha constrained).
* Running MCMC diagnostics for all the chains under each analysis.
* General discussion.

## Software

Before you start this tutorial, please make sure you have the following software installed on your PCs/HPCs:

* **`PAML`**: we used `PAML` v4.9i given that this was the most stable version when this project first started. Subsequent analyses until the end of this project therefore used the same version for consistency. If you have the latest PAML version, 4.10.7, available from the [`PAML` GitHub repository](https://github.com/abacus-gene/paml), you can run the programs aforementioned, but perhaps some values may change despite using the same seed numbers if there have been changes in the code from one version to another. Please see [the instructions on the PAML GitHub repository](https://github.com/abacus-gene/paml/wiki/Installation) to compile the version you are to use.

* **`R`** and **`RStudio`**: please download [R](https://cran.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/) so that you can run the various in-house R scripts used throughout the analyses detailed in this repository. The packages you will need to load should work with R versions that are either newer than or equal to v4.1.2 (we used v4.1.2). If you are a Windows user, please make sure that you have the correct version of `RTools` installed, which will allow you to install packages from the source code if required. For instance, if you have R v4.1.2, then installing `RTools4.0` shall be fine. If you have another R version installed on your PC, please check whether you need to install `RTools 4.2` or `RTools 4.3`. For more information on which version you should download, [please go to the CRAN website by following this link and download the version you need](https://cran.r-project.org/bin/windows/Rtools/).

    Before you proceed, however, please make sure that you install the following packages too:

    ```R
    # Run from the R console in RStudio
    # Check that you have at least R v4.1.2
    version$version.string
    # Now, install the packages we will be using
    # Note that it may take a while if you have not 
    # installed all these software before
    install.packages( c('rstudioapi', 'ape', 'phytools', 'sn', 'stringr', 'rstan', 'colorBlindness'), dep = TRUE )
    ## NOTE: If you are a Windows user and see the message "Do you want to install from sources the 
    ## packages which need compilarion?", please make sure that you have installed the `RTools`
    ## aforementioned.
    ## The versions we used for each of the packages aforementioned are the following:
    ##   rstudioapi: v0.14
    ##   ape: v5.7.1
    ##   phytools: v1.5.1
    ##   sn: v2.1.1
    ##   stringr: v1.5.0
    ##   rstan: v2.21.7
    ##   colorBlindness: v0.1.9
    ```

* **`FigTree`**: you can use this graphical interface to display tree topologies with/without branch lengths and with/without additional labels. You can then decide what you want to be displayed by selecting the buttons and options that you require for that to happen. You can [download the latest pre-compiled binaries, `FigTree v1.4.4`, from the `FigTree` GitHub repository](https://github.com/rambaut/figtree/releases).

* **`Tracer`**: you can use this graphical interface to visually assess the MCMCs you have run during your analyses (e.g., chain efficiency, chain convergence, autocorrelation, etc.). You can [download the latest pre-compiled binaries, `Tracer v1.7.2`, from the `Tracer` GitHub repository](https://github.com/beast-dev/tracer/releases/tag/v1.7.2).

## Citations

The dataset has been analysed using various tutorials, phylogenetics tools, and `MCMCtree` add-ons that Sandra ([@sabifo4](https://github.com/sabifo4)) is working on, which have been adapted to analyse the dataset used in this study. If you use the scripts provided here or in the [`Tutorials_MCMCtree` GitHub repository](https://github.com/sabifo4/Tutorial_MCMCtree), please cite the following:

**Álvarez-Carretero., S (2024). sabifo4/Tutorial_MCMCtree: v1.0.0 (tutorialMCMCtree-prerelease). Zenodo. https://doi.org/10.5281/zenodo.11306642**

> **NOTE**: once our manuscript has been peer-reviewed, we shall update this `README.md` file with more details with regards to how to cite our study and this repository.

## Data analysis

If you have gone through the previous sections and have a clear understanding of the dataset we used, the workflow we followed and you shall follow to reproduce our analyses, and have installed the required software to do so... Then you are ready to go!

You can start by taking a look at the input dataset we used prior to inferring the divergence times [by following this link](00_data/README.md).

Happy timetree inference! :)

----

Contact:

* Maria Eleonora Rossi: [`@MEleonoraRossi`](https://github.com/MEleonoraRossi).
* Sandra Álvarez-Carretero: [`@sabifo4`](https://github.com/sabifo4/).
