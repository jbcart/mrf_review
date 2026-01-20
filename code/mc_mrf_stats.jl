#!/usr/bin/env -S julia --threads=8

if length(ARGS) != 4
	error("Incorrect number of arguments.

	Usage: ./mc_mrf_stats.jl N ALPHA CODING NSIM
	")
end

using DataFrames, Statistics, StatsBase
using CSV
using ProgressBars
using Base.Threads
include("mrf_functions.jl")

# Settings
n = tryparse(Int64, ARGS[1])
α = tryparse(Float64, ARGS[2])
coding = ARGS[3]
nsim = tryparse(Int64, ARGS[4])

psis = 0:0.02:1.2 
# psis	= 0:.1:.1
P = length(psis)

ug = generatelattice(n,1)

N = n^2
n_edges = mrf_pairwise_ss(ug, ones(N), "ising")


mkpath("output/")

eT2i = zeros(P)
vT2i = zeros(P)
eT2a = zeros(P)
vT2a = zeros(P)
eT1 = zeros(P)
vT1 = zeros(P)
eDC = zeros(P)
vDC = zeros(P)

for p ∈ ProgressBar(eachindex(psis))
	T1 = zeros(nsim)
	T2i = zeros(nsim)
	T2a = zeros(nsim)
	DC = zeros(nsim)

	@threads for i ∈ 1:nsim
		if psis[p] > (0.82 + .5*abs(α)) && coding == "ising" && α != 0.0
			maxsweeps = 4
			M = 500
			method = "exact"
		elseif  α == 0.0 && coding == "ising"
			method = "sw"
			M = 400 
			maxsweeps = 0
		else
			maxsweeps = 0
			M = 250
			method = "exact"
		end
		config = rmrf(ug, α, psis[p], M,
			coding=coding, method=method, maxsweeps=maxsweeps, suppressmessage=true)

		T2i[i] = mrf_pairwise_ss(ug, config, "ising")
		T2a[i] = mrf_pairwise_ss(ug, config, "autologistic")
		T1[i] = sum(config)
		DC[i] = (countmap(config) |> values |> maximum)
# 		if psis[p] ∈ [.4,.6,.8,.84,.88,.92,.96,1.06,1.2]
# 			df = DataFrame(T1=T1, T2i=T2i, T2a=T2a, DC=DC)
# 			name = "output/" * coding[1:5] * "_n" * string(n)  *
# 			 	"_p" * replace(string(psis[p]), "." => "_") *
# 				"_a" * replace(string(α), "." => "_") * ".csv"
# 			CSV.write(name, df)
		end
	end

	eT1[p] = mean(T1)
	vT1[p] = var(T1)
	eT2i[p] = mean(T2i)
	vT2i[p] = var(T2i)
	eT2a[p] = mean(T2a)
	vT2a[p] = var(T2a)
	eDC[p] = mean(DC)
	vDC[p] = var(DC)
end
df = DataFrame(psi=psis, eT1=eT1, vT1=vT1,
 	eT2i=eT2i, vT2i=vT2i, eT2a=eT2a, vT2a=vT2a,
	eDC=eDC, vDC=vDC, N=N, n_edges=n_edges)

name = "output/" * coding[1:5] * "_n" * string(n) * "_a" * replace(string(α), "." => "_") * ".csv"
CSV.write(name, df)

