#!/usr/bin/env zsh

# Simulation runs to generate prior predictive response functions for direct data MRF models with an
# external field given by covariates.
./mc_mrf_grad.jl 2000 64 ising 0.0 0.0 1.0 int_norm
./mc_mrf_grad.jl 2000 64 autologistic 0.0 0.0 1.0 int_norm
./mc_mrf_grad.jl 2000 64 centered-auto 0.0 0.0 1.0 int_norm


# Simulation for response functions for the Ising model with alpha = 0.0 
./mc_mrf_stats.jl 64 0.0 ising 2000
