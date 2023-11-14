The dynamics of many gene regulatory networks is
characterized by the occurrence of metastable phenotypes and stochastic
phenotype switches. The chemical master equation (CME) is the most accurate
description to model such stochastic dynamics, whereby the long-time dynamics
of the system is encoded in the spectral properties of the CME operator. Markov
State Models (MSMs) provide a general framework for analyzing and visualizing
stochastic multistability and state transitions based on these spectral properties.
This code here implements a domain decomposition approach that
approximates the CME by a stochastic rate matrix on a discretized state space
and projects the multistable dynamics to a lower dimensional MSM. To
approximate the CME, the state space is decomposed via a Voronoi tessellation
and transition probabilities are estimated by using adaptive sampling strategies. 
The robust Perron cluster analysis (PCCA+) is applied to construct the final MSM.
Measures for uncertainty quantification are incorporated. 
The algorithm is applied to two GRN models, one for the
genetic toggle switch and one describing macrophage polarization.

# General Information

The code can be used to construct a Markov state model for two example systems:

model 1: macrophage polarization (https://doi.org/10.1016/j.jtbi.2023.111634) 

model 2: bistable toggle switch

as described in the paper

"Efficient construction of Markov state models for gene regulatory networks by domain decomposition"
by Maryam Yousefian, Anna-Simone Frank, Marcus Weber and Susanna Röblitz, arXiv 2023

Contributer: Maryam Yousefian (maryam.yousefian@uib.no), Susanna Röblitz (susanna.roblitz@uib.no)

Maintainer: Susanna Röblitz

# Included code files  and their description

The main file is `MSM.m`. It has been used to generate the files in the Results/ folder.
It adds the pathes PCCA/ and Dirichlet/.

The folder PCCA/ contains all necessary files to run the robust Perron cluster analysis.

The folder Dirichlet/ contains the two functions `dirichlet_fit_Newton.m` and `dirichlet_sample.m` from the 
fastfit toolbox (https://github.com/tminka/fastfit) as well as the two functions `digamma.m` and
`trigamma.m` from the lightspeed toolbox (https://github.com/tminka/lightspeed). These functions are needed
to estimate the parameters of the Dirichlet distribution and to sample from it. 

# Code outputs
Upon running `MSM.m`, two output mat-files are produced in the current folder:
One mat-file contains an object called `voronoi_table`, which contains all information about the discretization,
and an object called `Xh`, which contains all horizontal sampling points.
The second mat-file contains the matrix `alphaMat`, which has the estimated Dirichlet parameter vectors
alpha in its rows, and a matrix `Pmle`, which is the maximum likelihood estimate of the transition matrix.

# How to run the code
To reproduce the figures from the publication, run either `modelA_results.m` or `modelB_results.m`.
To construct a new MSM, modify the algorithmic parameters in `toggle_exp1.m` or `macrophage_exp3.m` and run `MSM.m`. 
