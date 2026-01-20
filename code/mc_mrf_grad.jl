#!/usr/bin/env -S julia --threads=8

if !(length(ARGS) ∈ [7])
	error("Incorrect number of arguments.

	Usage: ./mc_mrf_stats.jl NSIM N CODING ADJ BETA SD SCHEME
	")
end

using Distributions
using DataFrames, Statistics, StatsBase
using CSV
using ProgressBars
using Base.Threads
include("mrf_functions.jl")

# Settings
nsim = tryparse(Int64, ARGS[1])
n = tryparse(Int64, ARGS[2])
coding = ARGS[3]
adj = tryparse(Float64, ARGS[4])
β = tryparse(Float64, ARGS[5])
sd = tryparse(Float64, ARGS[6])
scheme = ARGS[7]
if sd != 0.0
  marginalize = true
else
  marginalize = false
end

if occursin("norm", scheme)
  prior_beta = Normal(β, sd)
elseif occursin("cauchy", scheme)
  prior_beta = Cauchy(β, sd)
elseif occursin("t7", scheme)
  prior_beta = β + TDist(7)*sd
end


method = "exact"

psi = 0:0.03:1.7
# psi	= 0:.1:.1
P = length(psi)
M = 250

ug = generatelattice(n,1)


N = n^2
n_edges = mrf_pairwise_ss(ug, ones(N), "ising")

df = DataFrame(x=repeat(1:n,inner=n), y=repeat(1:n,outer=n))
x_y = (df.x + df.y .- (n + 1)) ./ (n-1) .+ adj

α = zeros(N,2)

if !marginalize
	α[:,2] = x_y * β
	e_i = exp.(α[:,2]) ./ (1 .+ exp.(α[:,2]))
	config_est = Int.(e_i .> 0.5)
end

mkpath("output/")
name = "output/grad_" * coding[1:5] *
 	"_n" * string(n) *
	"_adj" * replace(string(adj), "." => "_") *
	"_b" * replace(string(β), "." => "_") *
	"_sd" * replace(string(sd), "." => "_") *
	scheme * ".csv"

eT2i = zeros(P)
vT2i = zeros(P)
uT2i = zeros(P)
lT2i = zeros(P)
eT2a = zeros(P)
vT2a = zeros(P)
uT2a = zeros(P)
lT2a = zeros(P)
eT1 = zeros(P)
vT1 = zeros(P)
uT1 = zeros(P)
lT1 = zeros(P)
eDC = zeros(P)
vDC = zeros(P)
uDC = zeros(P)
lDC = zeros(P)
eMAE = zeros(P)
uMAE = zeros(P)
lMAE = zeros(P)
eMR = zeros(P)
uMR = zeros(P)
lMR = zeros(P)


for p ∈ ProgressBar(eachindex(psi))
	T1 = zeros(nsim)
	T2i = zeros(nsim)
	T2a = zeros(nsim)
	DC = zeros(nsim)
	MAE = zeros(nsim)
	MR = zeros(nsim)
	config = zeros(N,nsim)
	if marginalize
		alpha = zeros(N,2,nsim)
	  e_imat = zeros(N,nsim)
		config_estmat = zeros(N,nsim)
	end

	@threads for i ∈ 1:nsim
		if psi[p] > 0.82 && coding ∈ ["centered-auto"]
			M = 500
			maxsweeps = 4
			method="exact"
		elseif psi[p] > 0.84 && coding ∈ ["ising"]
			method="sw-gibbs"
			M = 400
			maxsweeps = 0
		else
			M = 250 
			maxsweeps = 0
			method="exact"
		end
		if marginalize
			if occursin("int", scheme)
				alpha[:,2,i] .= x_y .* rand(prior_beta) .+ rand(prior_beta)
			else
				alpha[:,2,i] .= x_y .* rand(prior_beta)
			end
			e_imat[:,i] .= exp.(alpha[:,2,i]) ./ (1 .+ exp.(alpha[:,2,i]))
			config_estmat[:,i] .= Int.(e_imat[:,i] .> 0.5)

			config[:,i] = rmrf(ug, alpha[:,:,i], psi[p], M,
       	coding=coding, method=method, maxsweeps=maxsweeps, suppressmessage=true)
			
			MAE[i] = mean(abs.(config[:,i] .- e_imat[:,i]))
			MR[i] = mean(abs.(config[:,i] .- config_estmat[:,i]))
		else
			config[:,i] = rmrf(ug, α, psi[p], M,
      	coding=coding, method=method, maxsweeps=maxsweeps, suppressmessage=true)
			MAE[i] = mean(abs.(config[:,i] .- e_i))
			MR[i] = mean(abs.(config[:,i] .- config_est))
		end 

		T2i[i] = mrf_pairwise_ss(ug, config[:,i], "ising")
		T2a[i] = mrf_pairwise_ss(ug, config[:,i], "autologistic")
		T1[i] = sum(config[:,i])
		DC[i] = (countmap(config[:,i]) |> values |> maximum)
	end

	eT1[p] = mean(T1)
	vT1[p] = var(T1)
	uT1[p] = quantile(T1, 0.95)
	lT1[p] = quantile(T1,0.05)
	eT2i[p] = mean(T2i)
	vT2i[p] = var(T2i)
	uT2i[p] = quantile(T2i, 0.95)
	lT2i[p] = quantile(T2i, 0.05)
	eT2a[p] = mean(T2a)
	vT2a[p] = var(T2a)
	uT2a[p] = quantile(T2a, 0.95)
	lT2a[p] = quantile(T2a, 0.05)
	eDC[p] = mean(DC)
	vDC[p] = var(DC)
	uDC[p] = quantile(DC, 0.95)
	lDC[p] = quantile(DC, 0.05)
	eMAE[p] = mean(MAE)
	uMAE[p] = quantile(MAE, 0.95)
	lMAE[p] = quantile(MAE, 0.05)
	eMR[p] = mean(MR)
	uMR[p] = quantile(MR, 0.95)
	lMR[p] = quantile(MR, 0.05)
end

df = DataFrame(psi=psi, eT1=eT1, vT1=vT1, uT1=uT1, lT1=lT1,
 	eT2i=eT2i, vT2i=vT2i, uT2i=uT2i, lT2i=lT2i, eT2a=eT2a, vT2a=vT2a,
	uT2a=uT2a, lT2a=lT2a, eDC=eDC, vDC=vDC, uDC=uDC, lDC=lDC, eMAE=eMAE, 
	uMAE=uMAE, lMAE=lMAE, eMR=eMR, uMR=uMR, lMR=lMR,
	N=N, n_edges=n_edges)

CSV.write(name, df)


