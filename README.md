# mrf_review
Supplementary material for the paper "Markov Random Fields: Statistical Model Specification,
Phase Transition, and Response Functions."

## Response functions with a constant external field
### (Section 4.1)
The following executable file can be run to obtain Monte Carlo estimates of a set of response
functions for $\psi \in \{0.0, 0.02, 0.04, \dots, 1.2\}$ for the three binary formulations of an MRF
with a constant external field (i.e., $\alpha = C$, where $C$ is a constant).

./mc_mrf_stats.jl N ALPHA CODING NSIM

N       - The number of rows in an NxN first-order lattice.
ALPHA   - The value of the external field parameter $\alpha$.
CODING  - The formulation of the MRF, can be "autologistic", "ising", or "ising-physics".
NSIM    - The number of Monte Carlo simulations for each value of $\psi$. 

The command to generate the response functions in the paper is:

./mc_mrf_stats.jl 64 0.0 ising 2000

\**Note: the above executable file specifies the number of threads for julia to use in the shebang line*

## Prior predictive response functions for direct data models with a covariate external field.
### (Section 5.1.2)
The following executable file can be run to obtain Monte Carlo estimates of a set of response
functons for $\psi \in \{0.0, 0.03, 0.03, \dots, 1.7\}$ for the autologistic, centered-autologistic 
and Ising model with an external field of covariates $x_i=(r(i)+c(i)-n_r-1)/(n_r-1)$,
where $r(i)$ returns the row of areal unit $i$, $c(i)$ returns the column, and $n_r=64$, the number
of rows/columns. The external field is given by $\alpha_i = \beta_0 + (x_i+a)\beta_1$ and $\beta_0$, 
$\beta_1$ are distributed normal, cauchy, or $t_7$ with scale parameter $\sigma$ and location $\mu$. 

./mc_mrf_stats.jl NSIM N CODING ADJ BETA SD SCHEME

| Parameter | Description
| --- | ---
| NSIM    | The number of Monte Carlo simulations for each value of $\psi$.
| N       | The number of rows in an NxN first-order lattice.
| CODING  | The formulation of the MRF, can be "autologistic", "ising", or "centered-auto".
| ADJ     | An adjustment to the gradient of covariates (the $a$ in $x_i+a$ above).
| BETA    | The location parameter for the prior on $\beta_0$ and $\beta_1$.
| SD      | The scale parameter for the prior on $\beta_0$ and $\beta_1$.
| SCHEME  | A string tag to indicate the prior on $\beta$  (cauchy, t7 or, norm) and whether to include an intercept (int) in the model for the external field. For example "int_norm" will include an intercept $\beta_0$ (e.g., $\alpha_i = \beta_0 + (x_i+a)\beta_1$) and $\beta \sim \mathrm{norm}(\text{BETA}, \text{SD})$.

The commands to generate the prior predictive response functions in the paper are:

./mc_mrf_grad.jl 2000 64 ising 0.0 0.0 1.0 int_norm
./mc_mrf_grad.jl 2000 64 autologistic 0.0 0.0 1.0 int_norm
./mc_mrf_grad.jl 2000 64 centered-auto 0.0 0.0 1.0 int_norm

\**Note: the above executable files specify the number of threads for julia to use in the shebang line*

## Run Simulations
All of the above simulations can be run with:

./run_simulations.zsh

Remember to specify the number of threads in both mc_mrf_stats.jl and mc_mrf_grad.jl before running.

## Plots
Code to recreate all figures in the paper is in mrf_review_plots.r
