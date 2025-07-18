```text
TITLE: Repository of "Media Bias and Polarization through the Lens of a Markov Switching Latent Space Network Model"

AUTHORS:      Roberto Casarin, Antonio Peruzzi, Mark FJ Steel

AVAILABLE AT:    The Annals of Applied Statistics

PLEASE CITE AS:   (forthcoming)

DATE:       June 2025


Tested on R version 4.2.2 (2022-10-31) -- "Innocent and Trusting"

-------------------------------------------------------------------------------
The following Repository contains the files (scripts and data) used to reproduce 
the results of the paper "Media Bias and Polarization through the Lens of a Markov 
Switching Latent Space Network Model".
-------------------------------------------------------------------------------

        %%%%%%%%%%%%%% PRELIMINARIES  %%%%%%%%%%%%%%

-------------------------------------------------------------------------------

Our MCMC algorithm is entirely implemented in C++, enabling faster execution 
speed compared to interpreted languages like R or Python. However, we still rely on R 
for data manipulation and plotting. The smooth integration of the two languages has 
been made possible through the utilization of the Rcpp package, which offers a convenient 
interface for invoking C++ scripts within R.

Before running the following scripts, make sure that your version of R is updated (at least v. 2020-10-10)
and to run the script Preliminary_01_InstallPackages.R inside the folder "01-Preliminaries" to
install all the required packages. Make also sure to change the working directory to your path whenever 
It is clearly stated:


########CHANGE YOUR PATH ###########
setwd("~/Desktop/Repository/")
#####################################

-------------------------------------------------------------------------------

        %%%%%%%%%%%%%%  R SCRIPTS  %%%%%%%%%%%%%%

-------------------------------------------------------------------------------


The following R script files are used to estimate the Bayesian Markov-Switching Latent-Space
network model (hereafter MS-LS) on the datasets studied in the main paper.
(Section 3, Section 4, and Supplement):

* Dynamic_01_Results_[Country].R
Estimates the MS-LS model on the dynamic network dataset.
Running time (*) > 20 hrs 

* Static_01_Results_[Country].R
Estimates the MS-LS model on the static network dataset.
Running time ~45 mins

* Simulation_02_results.R
Generates a synthetic dataset and estimates the MS-LS model.
Running time ~13 mins



The following R scripts files allow us to plot the Figures of a MS-LS model with D = 1 and K = 2:

* Graph_01.R
Generates the introductory Graph reported in Figure 1 and Figure 2

* Properties_01.R
Generates the contour plots reported in Figure I.1

* Simulation_02_results.R
Generates the figures of the synthetic data analysis (Figure 5, Figure D.1, Figure D.2, Figure D.3, Table D.1)

* Static_02_Plots.R
Generates the figures of the static data analysis (Figure 7, Figure I.1)

* Dynamic_02_Plots.R
Generates the figures of the dynamic data analysis for an MS-LS with d = 1 and K = 2
(Figure I.1, Figure I.2, Figure I.3, Figure I.4).

The  R scripts for running an MS-LS model with D > 1 and K > 2
are contained in the Folder named "Extended".

[Country] = DE,FR,IT,SP (*) with MacBook Air 2022 M2 

-------------------------------------------------------------------------------

        %%%%%%%%%%%%%%  C++ Scripts  %%%%%%%%%%%%%%

-------------------------------------------------------------------------------

The C++ script necessary for running the Bayesian MS-LS model is reported here below along
with a brief description of the main function and its signature.C++ is integrated into R via Rcpp.

* MS_LS_FE.cpp
The script contains the MCMC function to estimate the Bayesian MS-LS network model for a dynamic network.
The dynamic network is expected to have N nodes for each time t=1,...,Time and count weighted edges.


INPUT

The MCMC function expects the use of the columns of 3 data frames as main input: "EL_x", "EL_princ", "DBplane".

* "EL_x" is an edge list of the off-diagonal elements of a network through time with columns:
  -"i": node i
  -"j": node j
  -"t": time t
  -"w": countable weight 
  "one": value 1

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~ example with N = 3 and T = 2 ~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


----------------------------
 "i", "j", "t", " w", "one"
----------------------------
 1  2   1  w121  1
 1  3   1  w131  1
 2  1   1  w211  1
 2  3   1  w231  1
 3  1   1  w311  1
 3  2   1  w321  1
 1  2   2  w122  1
 1  3   2  w132  1
 2  1   2  w212  1
 2  3   2  w232  1
 3  1   2  w312  1
 3  2   2  w322  1
----------------------------
  
* "EL_princ" is an edge list of the lower-triangle elements (i > j) of a network through time with columns:
  -"i": node i
  -"j": node j
  -"t": time t
  -"w": countable weight 
  "one": value 1

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~ example with N = 3 and T = 2 ~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


----------------------------
 "i", "j", "t", " w", "one"
----------------------------
 1  2   1  w121  1
 1  3   1  w131  1
 2  3   1  w231  1
 1  2   2  w122  1
 1  3   2  w132  1
 2  3   2  w232  1
----------------------------


* "DBplane" is a dataset of nodes' features through time with columns:
  -"i": node i
  -"t": time t
  -"leaning": leaning index 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~ example with N = 3 and T = 2 ~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


---------------------
 "i", "t", "leaning",
---------------------
 1  1   l11
 2  1   l21
 3  1   l31
 1  2   l12
 2  2   l22
 3  2   l23
---------------------


FUNCTION SIGNATURE
Here below is a description of all the arguments of the function MCMC

result = MCMC(princ_w = EL_princ$w, ..................#edge weights from EL_princ
       princ_ones = EL_princ$ones, ............#vector of ones from from EL_princ
       x_w = EL_x$w, ..........................#edge weights from EL_x
       x_ones = EL_x$ones, ...................#vector of ones from EL_x
       leaning = DBplane$leaning, .............#text analysis leaning from DBplane
       lam_ad_beta = lam_ad_beta, .............#adaptive RW-MH lambda parameter - individual effects
       mu_mat_beta = mu_mat_beta, .............#adaptive RW-MH mu - individual effects
       Sigma_ad_beta = Sigma_ad_beta, .........#adaptive RW-MH sigma - individual effects
       lam_ad_za = lam_ad_za, .................#adaptive RW-MH lambda - latent coordinate - state a (vector of length N)
       mu_mat_za = mu_mat_za, .................#adaptive RW-MH mu - latent coordinate - state a (matrix of dimension Nx2)
       Sigma_ad_za = Sigma_ad_za, .............#adaptive RW-MH sigma - latent coordinate - state a (list of N 2x2 matrices)
       lam_ad_zb = lam_ad_zb, .................#adaptive RW-MH lambda - latent coordinate - state b (see above)
       mu_mat_zb = mu_mat_zb, .................#adaptive RW-MH mu - latent coordinate - state b (see above)
       Sigma_ad_zb = Sigma_ad_zb, .............#adaptive RW-MH sigma - latent coordinate - state b (see above)
       beta = rnorm(Npages, 0, 0.01), ........#starting value - individual effect (vector of length N)
       xi_state1 = xi[,1], ....................#starting value - state 1 (boolean vector of length Times)
       zi_a_1 = rep(0, Npages), ...............#starting value - latent coordinate - state a (vector of length N)
       zi_b_1 = rep(0, Npages), ...............#starting value - latent coordinate - state b (vector of length N)
       mu_beta =0, ............................#prior mean - individual effect
       sigma_beta = 15, .......................#prior sd - individual effect
       mu_a = 0, .............................#prior mean - latent coordinate - state a
       sigma_a = 10, ..........................#starting sd - latent coordinate - state a
       mu_b = 0, .............................#prior mean - latent coordinate - state b
       sigma_b = 10, ..........................#starting sd - latent coordinate - state b
       s_a = 0.1,..............................#prior shape parameter of sigma (gamma(s_a, s_b))
       s_b = 0.1, .............................#prior scale parameter of sigma (gamma(s_a, s_b))
       phi = 50, ..............................#starting value phi
       gamma_0 = 0, ...........................#starting value - gamma_0
       gamma_1 = 0, ...........................#starting value - gamma_1
       a_phi = a_phi, .........................#prior shape parameter of phi (gamma(a_phi, b_phi))
       b_phi = b_phi, .........................#prior scale parameter of phi (gamma(a_phi, b_phi))
       a_gamma_0 = 0, .........................#prior mean parameter of gamma_0
       b_gamma_0 = 15, ........................#prior sd parameter of gamma_0
       a_gamma_1 = 0, .........................#prior mean parameter of gamma_1
       b_gamma_1 = 15, ........................#prior sd parameter of gamma_1
       omega_lower_a = 2, .....................#starting value - 1st parameter dirichlet
       omega_lower_b = 2, .....................#starting value - 2nd parameter dirichlet
       P = P, .................................#2x2 transition Matrix P
       N = Npages, ............................# Number of nodes
       Time = 100, ............................#Number of times
       x_i = EL_x$i, ..........................#index vector of nodes in EL_x
       DBplane_i = DBplane$i, ................#index vector of nodes in DBplane
       prop_sd_gamma = 0.01, ..................#proposal RWMH for the gamma parameters
       prop_sd_phi = 35, ......................#proposal RWMH for the phi parameter
       acc_beta = 0.25, .......................#target acceptance Adaptive RW-MH - individual effects
       acc_zeta_a = 0.25, .....................#target acceptance Adaptive RW-MH - latent coordinates - state a
       acc_zeta_b = 0.25, .................... #target acceptance Adaptive RW-MH - latent coordinates - state b
       pivot = 3, .............................#pivot node (news outlet) to solve identification
       sign= -1, ..............................#sign of the news outlet -1 left or 1 right
       interp_eq = interp_eq, .................#0-1 option for the use of the interpretation equation
       ms_eq = ms_eq, .........................#0-1 option for the use of the markov switching
       rg_eq = rg_eq,..........................#0-1-2 option random graph
       Iterations = Iterations.................#number of iterations
       )

NOTE on the use of options:

- interp_eq = 0 implies that the MS-LS model is computed disregarding equation 3 (see the manuscript)
- ms_eq = 0 implies that the LS model is computed diregarding the Markov-Switching component
- rg_eq = 1 implies that a simple random graph model y_{ijt} ~ Pois(exp(alpha)) is run
- rg_eq = 2 implies that a random graph model y_{ijt} ~ Pois(exp(alpha_i + alpha_j )) is run
- the combination interp_eq = 1, ms_eq = 1, rg_eq = 0 returns the MS-LS model as in the manuscript

OUTPUT

The main output is a list object called "result" containing the following elements:

- result[[1]] the matrix beta_it containing the iteration draws for the individual effects
- result[[2]] the matrix phi_gamma_it containing the iteration draws for the parameters gamma_0, gamma_1, phi
- result[[3]] the matrix zi_a_it containing the iteration draws for the parameters zi_a, the latent coordinates in state a
- result[[4]] the matrix zi_b_it containing the iteration draws for the parameters zi_b, the latent coordinates in state b
- result[[5]] the matrix P_ite containing the iteration draws for the parameters p1 and p2
- result[[6]] the matrix x_it containing the iteration draws for the latent states 
- result[[7]] the matrix HETA containing the likelihood in state a, in state b and given the current x_it
- result[[8]] the matrix Sigma_z_ite containing the iteration draws for the variances of the latent states 

-------------------------------------------------------------------------------

        %%%%%%%%%%%%%%  Other R SCRIPTS  %%%%%%%%%%%%%%

-------------------------------------------------------------------------------


* 01_TextAnalysis_IT.R
We provide this script as an example of how to create the "slant index":
The code implements a time-varying cosine similarity between the language used by political parties
And the language used by news outlets on a daily basis. The output is the DBplane dataset for Italy.
Running time ~45 mins.

The code is also used to reproduce Figure F.1 and F.2. 

* Properties_02.R
Generates the contour plots reported in Figure B.2
Running time ~10 mins.



-------------------------------------------------------------------------------

        %%%%%%%%%%%%%%  Other C++ SCRIPTS  %%%%%%%%%%%%%%

-------------------------------------------------------------------------------

The following script can be used to run a static version of the MS-LS model:

* Rcpp_rf_[country].cpp
The script contains a function, namely MCMC, suitable for static analysis (similar to the main c++ script).

* Predictive.cpp
The script contains a set of functions to generate Table 2, Table 4 and Figure I.2


-------------------------------------------------------------------------------

         %%%%%%%%%%%%%%  DATA  %%%%%%%%%%%%%%

-------------------------------------------------------------------------------

The following .RData files include the data used in the applications (Section 4):

* SimulationEnv_FE.RData
Contains the synthetic dataset used in Section 3.3 

* Data_Env_single_[country].RData
Contains the static dataset used in Section 4.2

 + we refer to the datasets EL_x and EL_princ as "network dataset" in the manuscript
 + we refer to the dataset DBplane as "media slant index" in the manuscript

* DataEnv_[country]_all.RData
Contains the dynamic dataset used in Section 4.3

 + we refer to the datasets EL_x and EL_princ as "network dataset" in the manuscript
 + we refer to the dataset DBplane as "media slant index" in the manuscript

DATA CONSTRUCTION 

* We created the "network dataset" from the Facebook dataset provided by:

Schmidt, A. L., F. Zollo, A. Scala, and W. Quattrociocchi (2018). Polarization Rank: A
Study on European News Consumption on Facebook. arXiv preprint arXiv:1805.08030.

Our time-varying networks are briefly described here below.
Nodes represent news outlets and the weight associated with each edge represents the number of "commenters" in common between
Node i and node j at time t for t = 1:Time

 + Germany Nodes = 47 Time = 729 days; from "2015-01-01" to "2016-12-31"
 + France Nodes = 62 Time = 729 days;
 + Italy  Nodes = 45 Time = 729 days;
 + Spain  Nodes = 43 Time = 729 days;

* We created the "slant index" from CrowdTangle Facebook data.
An illustrative sample is available following the path Data/TextAnalysis-SlantIndex/Italy

-------------------------------------------------------------------------------

        %%%%%%%%%%%%%%  EXTENDED  %%%%%%%%%%%%%%

-------------------------------------------------------------------------------

The "Extended" folder contains C++ and R scripts designed to estimate the MS-LS model
in a more flexible setting, allowing for an arbitrary number of Markov states \( K \) and latent dimensions \( D \).

This version is suitable for more complex applications, such as higher-dimensional embeddings of media outlets
or the presence of more than two polarization regimes. It generalizes the basic version described above.

The main C++ script for the extended model is:

* Extended_MS_LS_FE.cpp  
  Contains the MCMC function for estimating the Bayesian MS-LS model with flexible dimension \( D \) and number of regimes \( K \).  
  Integrated into R via **Rcpp**.

**FUNCTION SIGNATURE**  
Below is the complete signature of the MCMC function implemented in this script:

result = MCMC(princ_w = EL_princ$w, ........................# edge weights from the network triangle dataset
              princ_ones = EL_princ$ones,...................# vector of ones from the triangle dataset
              x_w = EL_x$w,.................................# edge weights from the network off-diagonal dataset
              x_ones = EL_x$ones,,..........................# vector of ones from the off-diagonal dataset
              leaning = DBplane$leaning,................... # text analysis leaning (media slant index)
              lam_ad_beta = lam_ad_beta,....................# adaptive RW-MH lambda - individual effects
              mu_mat_beta = mu_mat_beta,....................# adaptive RW-MH mu - individual effects
              Sigma_ad_beta = Sigma_ad_beta,................# adaptive RW-MH sigma - individual effects
              lam_ad_z = lam_ad_z,..........................# adaptive RW-MH lambda - latent positions
              mu_mat_z = mu_mat_z,..........................# adaptive RW-MH mu - latent positions
              Sigma_ad_z = Sigma_ad_z,......................# adaptive RW-MH sigma - latent positions
              beta = rnorm(Npages, 0 , 0.01),...............# starting value - individual effects
              xi_state = as.matrix(xi[, 1:K]),..............# starting regime assignments
              zi_in = zi,...................................# starting value - latent positions (array N x D)
              mu_beta = 0, .................................# prior mean - individual effect
              sigma_beta = 15,..............................# prior sd - individual effect
              mu_z = rep(0, D),.............................# prior mean - latent coordinates
              sigma_z = rep(10, K),........................ # prior sd - latent coordinates for each regime
              s_a = 0.01, s_b = 0.01,.......................# prior shape and rate parameters for latent coordinate variance
              phi = 50,.....................................# starting value - phi
              gamma_0 = 0, gamma_1 = 1,.....................# starting values - gamma parameters
              a_phi = a_phi, b_phi = b_phi,.................# prior hyperparameters for phi
              a_gamma_0 = 0, b_gamma_0 = 15,................# prior hyperparameters for gamma_0
              a_gamma_1 = 0, b_gamma_1 = 15,................# prior hyperparameters for gamma_1
              omega = rep(1/K, K),..........................# initial state probabilities (Dirichlet)
              P = P,........................................# K x K transition matrix
              N = Npages,...................................# number of nodes
              D = D,........................................# latent space dimension
              K = K,........................................# number of Markov regimes
              Time = Time,..................................# number of time periods
              x_i = EL_x$ith,...............................# index vector of nodes in EL_x
              DBplane_i = DBplane$i,........................# index vector of nodes in DBplane
              prop_sd_gamma = 0.01,.........................# RW-MH proposal sd for gamma parameters
              prop_sd_phi = 10,.............................# RW-MH proposal sd for phi
              acc_beta = 0.25,..............................# target acceptance - individual effects
              acc_zeta = 0.25...............................# target acceptance - latent positions
              pivot = which(pages_names == "xx"),...........# pivot node for identification
              sign = 1,.....................................# sign for identification constraint
              interp_eq = 1,................................# include interpretation equation (0/1)
              ms_eq = 1,....................................# include Markov-switching dynamics (0/1)
              rg_eq = 0,................................... # graph model choice (0/1/2)
              Iterations = Iterations,......................# number of MCMC iterations
)


```
