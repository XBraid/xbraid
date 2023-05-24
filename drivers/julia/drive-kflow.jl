using ForwardDiff, DiffResults, LinearAlgebra, PreallocationTools
using Interpolations, FFTW
using MPI, Base.Threads
using BenchmarkTools, Plots, LaTeXStrings
using Statistics: mean
using Random: seed!
using LinearAlgebra: norm

include("../../braid/braid.jl/XBraid.jl")
using .XBraid

MPI.Init()
comm = MPI.COMM_WORLD

Float = Float64

"""
struct containing everything needed internally for the simulation
"""
struct KFlowApp
	cf::Int
	useTheta::Bool
	deltaRank::Int
	coords::Array{Float, 3}
	κ::Frequencies{Float}
	k_num::Int
	ℜ::Float
	solution::Vector{Array{Float, 3}}
	solTf::Vector{Array{Float, 3}}
	lyapunov_vecs::Vector{Vector{Array{Float, 3}}}
	lyapunov_exps::Vector{Vector{Float}}
	times::Vector{Int}
	x_d::DiffCache{Array{Float, 3}}
	u_d::DiffCache{Array{Float, 3}}
	ϕ_d::DiffCache{Array{Complex{Float}, 2}}
	P̂::FFTW.FFTWPlan
end

# default constructor
function KFlowApp(cf::Integer, useTheta::Bool, deltaRank::Integer, k_num::Integer, ℜ::Real)
	coords = zeros(Float, nₓ, nₓ, 2)
	x_range = range(0, lengthScale - Δx, nₓ)
	@views @inbounds for (i, x) ∈ enumerate(x_range), (j, y) ∈ enumerate(x_range)
		coords[i, j, 1] = x
		coords[i, j, 2] = y
	end
	κ = 2π / lengthScale .* fftfreq(nₓ, nₓ)
	# preallocations
	x_d = zeros(Float, nₓ, nₓ, 2)
	u_d = zeros(Float, nₓ, nₓ, 2)
	ϕ_d = zeros(Complex{Float}, nₓ, nₓ)
	KFlowApp(cf, useTheta, deltaRank, coords, κ, k_num, ℜ, [], [], [], [], [], DiffCache(x_d), DiffCache(u_d), DiffCache(ϕ_d), plan_fft(coords, [1, 2]))
end

oneball(n) = [(n + 1 - i, i) for i ∈ 1:n]
wavenumbers(shell) = reduce(hcat, [oneball(i) for i ∈ 1:shell])
is_in_shell(i, shell) = i <= trunc(Int, shell * (shell + 1) / 2)

"""
compute ith wavenumber for a 2D scalar field
"""
function get_wavenumber2D(i)
	shell = 1
	for _ ∈ 1:i
		is_in_shell(i, shell) && break
		shell += 1
	end
	shell_start = trunc(Int, shell * (shell - 1) / 2)
	rel_i = i - shell_start
	return oneball(shell)[rel_i]
end

"""
compute the ith Fourier mode of a 1D scalar field
"""
function fourierMode1D(x::Real, k::Integer)
	k % 2 == 1 && return cos(trunc(k / 2) * 2π / lengthScale * x)
	k % 2 == 0 && return sin(trunc((k + 1) / 2) * 2π / lengthScale * x)
end

"""
compute the ith Fourier mode of a 2D scalar field
"""
function fourierMode2D!(app::KFlowApp, u_field, i::Integer)
	kx, ky = get_wavenumber2D(i)
	@views x, y = app.coords[:, :, 1], app.coords[:, :, 2]
	@. u_field = fourierMode1D(x, kx) * fourierMode1D(y, ky)
end

"""
compute the ith Fourier mode of a 2D vector field
"""
function fourierMode2DVec!(app::KFlowApp, u, i::Integer)
	@views ux, uy = u[:, :, 1], u[:, :, 2]
	kx, ky = get_wavenumber2D(i)
	fourierMode2D!(app, ux, kx)
	fourierMode2D!(app, uy, ky)
end

"""
set the initial condition to a Taylor-Green vortex
"""
function TaylorGreen!(app::KFlowApp, u, k = 1)
	coords = app.coords
	kx = trunc(Int, k / 2) + k % 2
	ky = trunc(Int, k / 2) + 1
	@views u_x, u_y = u[:, :, 1], u[:, :, 2]
	@views x, y = coords[:, :, 1], coords[:, :, 2]
	@. u_x = sin(kx * x) * cos(ky * y)
	@. u_y = -cos(kx * x) * sin(ky * y)
	return
end

"""
compute the forcing term
"""
function kolmogorovForce!(app::KFlowApp, u, Δt; μ = 32.)
	k = app.k_num
	# Fₖ(x, y) = (sin(ky), 0)
	@views @. u[:, :, 1] += μ * Δt * sin(k * app.coords[:, :, 2])
	return u
end

"""
compute self advection of u
"""
function advect_semi_lagrangian!(app::KFlowApp, u::AbstractArray, v::AbstractArray, Δt::Real)
	# get cached arrays (either real or dual, depending on typeof(u))
	coords_new = get_tmp(app.x_d, u)
	u_new = get_tmp(app.u_d, u)

	# (1) trace each grid point backwards along u
	@. coords_new = app.coords - Δt * v

	# (2) interpolate the velocity field at the new grid points
	x_range = range(0.0, lengthScale - Δx, nₓ)
	for i in 1:2
		itp = interpolate(u[:, :, i], BSpline(Linear(Periodic())))
		sitp = scale(itp, x_range, x_range)
		extp = extrapolate(sitp, Periodic())

		u_new[:, :, i] .= @views extp.(coords_new[:, :, 1], coords_new[:, :, 2])
	end
	u .= u_new
	return u
end

"""
Solve Poisson problem for incompressible velocity field
Δϕ = ∇⋅u
u = u - ∇ϕ
"""
function project_incompressible!(app::KFlowApp, u::AbstractArray, Δt::Real)
	κ = get_tmp(app.x_d, u)
	ϕ = get_tmp(app.ϕ_d, u)
	P̂ = app.P̂
	ℜ = app.ℜ
	for (i, κ_x) ∈ enumerate(app.κ), (j, κ_y) ∈ enumerate(app.κ)
		κ[i, j, 1] = κ_x
		κ[i, j, 2] = κ_y
	end
	û = P̂ * u
	@views κ_x, κ_y = κ[:, :, 1], κ[:, :, 2]
	@views û_x, û_y = û[:, :, 1], û[:, :, 2]
	# diffusion
	@. û_x *= exp(-Δt / ℜ * (κ_x^2 + κ_y^2))
	@. û_y *= exp(-Δt / ℜ * (κ_x^2 + κ_y^2))

	# KS equation: uₜ + u(∇u) = - Δu - Δ²u
	# @. û_x *= exp(Δt*(κ_x^2 + κ_y^2 - (κ_x^2 + κ_y^2)^2))
	# @. û_y *= exp(Δt*(κ_x^2 + κ_y^2 - (κ_x^2 + κ_y^2)^2))

	# incompressible projection
	# ϕ = inv(Δ)∇⋅u
	@. ϕ = -(1.0im * κ_x * û_x + 1.0im * κ_y * û_y) / (κ_x^2 + κ_y^2)
	ϕ[1, 1] = 0.0 + 0.0im # pressure is unique up to a constant
	# u = u - ∇ϕ
	@. û_x -= 1.0im * κ_x * ϕ
	@. û_y -= 1.0im * κ_y * ϕ

	u .= real(P̂ \ û)
	return u
end

function extractPartials(u::AbstractArray{ForwardDiff.Dual{T, V, P}}) where {T, V, P}
	ps = zeros(V, size(u)..., P)
	for I ∈ CartesianIndices(u), j ∈ 1:P
		@inbounds ps[I, j] = ForwardDiff.partials(u[I], j)
	end
	return ps
end

function fillArrayOfDuals!(u::AbstractArray{ForwardDiff.Dual{T, V, P}}, vs::AbstractArray{V}, ps::AbstractArray{V}) where {T, V, P}
	checkbounds(ps, first(CartesianIndices(u)), 1:P) # make inbounds safe
	for I ∈ CartesianIndices(u)
		@inbounds u[I] = ForwardDiff.Dual{T}(vs[I], ntuple(j -> @inbounds(ps[I, j]), P))
	end
end

# this enables ForwardDiff through the FFT where it normally doesn't work
function project_incompressible!(app::KFlowApp, u::AbstractArray{ForwardDiff.Dual{T, V, P}}, Δt::Real) where {T, V, P}
	vs = ForwardDiff.value.(u)
	ps = extractPartials(u)
	project_incompressible!(app, vs, Δt)
	map(eachslice(ps, dims = 4)) do p
		project_incompressible!(app, p, Δt)
	end
	fillArrayOfDuals!(u, vs, ps)
	return u
end

function ∇_dot(app, u)
	κ = get_tmp(app.x_d, u)
	for (i, κ_x) ∈ enumerate(app.κ), (j, κ_y) ∈ enumerate(app.κ)
		κ[i, j, 1] = κ_x
		κ[i, j, 2] = κ_y
	end
	@views κx, κy = κ[:, :, 1], κ[:, :, 2]
	û = app.P̂ * u
	out = @. 1.0im * κx * û[:, :, 1] + 1.0im * κy * û[:, :, 2]
	return real(ifft(out))
end

function ∇X(app::KFlowApp, V)
	κ = get_tmp(app.x_d, V)
	for (i, κ_x) ∈ enumerate(app.κ), (j, κ_y) ∈ enumerate(app.κ)
		κ[i, j, 1] = κ_x
		κ[i, j, 2] = κ_y
	end
	f = get_tmp(app.ϕ_d, V)
	V̂ = app.P̂ * V
	@views V̂x, V̂y = V̂[:, :, 1], V̂[:, :, 2]
	@views κx, κy = κ[:, :, 1], κ[:, :, 2]
	out = @. 1.0im * κx * V̂y - 1.0im * κy * V̂x
	return real(ifft(out))
end

function spectralInterpolation(f, targetSize)
	size_f = size(f)[1] # assume square
	targetSize > size_f || return f
	npad = floor(Int, targetSize / 2 - size_f / 2)
	f̂ = fftshift(fft(f))
	C = eltype(f̂)
	f̂_interp = [
		zeros(C, npad, npad)   zeros(C, npad, size_f) zeros(C, npad, npad)
		zeros(C, size_f, npad) f̂           zeros(C, size_f, npad)
		zeros(C, npad, npad)   zeros(C, npad, size_f) zeros(C, npad, npad)
	]
	return real(ifft(ifftshift(f̂_interp)))
end

function base_step!(app::KFlowApp, u::AbstractArray, Δt::Float)
	kolmogorovForce!(app, u, Δt)
	advect_semi_lagrangian!(app, u, u, Δt)
	project_incompressible!(app, u, Δt)
	return u
end

function my_init(app::KFlowApp, t::Float)
	seed!(1)
	# u = zeros(Float, nₓ, nₓ, 2)
	u = randn(Float, nₓ, nₓ, 2)
	# TaylorGreen!(app, u)
	project_incompressible!(app, u, 1e-1app.ℜ)
	u[:, :, 1] .-= mean(u[:, :, 1])
	u[:, :, 2] .-= mean(u[:, :, 2])
	u ./= my_norm(app, u)
	return u
end

function my_basis_init(app::KFlowApp, t::Float, i::Integer)
	ψ = zeros(Float, nₓ, nₓ, 2)
	fourierMode2DVec!(app, ψ, i + 1)
	return ψ
end

function my_step!(
	app::KFlowApp, status::Ptr{Cvoid},
	u::AbstractArray, ustop::AbstractArray,
	tstart::Real, tstop::Real,
)
	Δt = tstop - tstart
	level = XBraid.status_GetLevel(status)
	iter = XBraid.status_GetIter(status)
	if !app.useTheta || level == 0
		u = base_step!(app, u, Δt)
	else
		# richardson based θ method
		m = app.cf ^ level
		p = 1
		θ = 2^p * (m^p - 1) / (m^p * (2^p - 1))
		# θ = 2
		u_sub = deepcopy(u)
		base_step!(app, u_sub, Δt / 2)
		base_step!(app, u_sub, Δt / 2)
		base_step!(app, u, Δt)
		@. u = θ * u_sub + (1 - θ) * u
	end
	return u
end

function my_step!(
	app::KFlowApp, status::Ptr{Cvoid},
	u::AbstractArray, ustop::AbstractArray,
	tstart::Real, tstop::Real,
	Ψ::Vector{Any}
)
	rank = length(Ψ)
	Ψ_new = reduce((a, b) -> cat(a, b, dims = 4), Ψ)
	perturb(r) = my_step!(app, status, u + r' * Ψ, ustop, tstart, tstop)

	result = DiffResults.DiffResult(u, Ψ_new)
	result = ForwardDiff.jacobian!(result, perturb, zeros(rank))
	for i in 1:rank
		Ψ[i] .= Ψ_new[:, :, :, i]
	end
	return
end

function my_sum!(app, a, x, b, y)
	@. y = a * x + b * y
end

function my_access(app::KFlowApp, status, u)
	XBraid.status_GetWrapperTest(status) && return
	ntime = XBraid.status_GetNTPoints(status)
	index = XBraid.status_GetTIndex(status)
	level, done = XBraid.status_GetLevel(status), XBraid.status_GetDone(status)

	if level == 0 && done && index % app.cf == 0
		push!(app.solution, deepcopy(u))
		push!(app.times, index)
		rank = XBraid.status_GetDeltaRank(status)
		if rank > 0
			exps = XBraid.status_GetLocalLyapExponents(status)
			vecs = XBraid.status_GetBasisVectors(status)
			push!(app.lyapunov_exps, exps)
			push!(app.lyapunov_vecs, deepcopy(vecs))
		end
	end

	if level == 0 && index == ntime-1
		# last time step
		push!(app.solTf, deepcopy(u))
	end
end

my_norm(app, u) = LinearAlgebra.normInf(u)
my_innerprod(app, u, v) = u ⋅ v

const lengthScale = 2π
const nₓ = 99
const Δx = lengthScale / nₓ

function test()
	my_app = KFlowApp(2, false, 3, 2, 1600)

	test_app = XBraid.BraidApp(
		my_app, comm, comm,
		my_step!, my_init,
		my_sum!, my_norm, my_access,
		my_basis_init, my_innerprod)

	XBraid.testInitAccess(test_app, 0.0)
	XBraid.testClone(test_app, 0.0)
	XBraid.testSpatialNorm(test_app, 0.0)
	XBraid.testBuf(test_app, 0.0)
	XBraid.testDelta(test_app, 0.0, 0.1, my_app.deltaRank)

	# curl = ∇X(my_app, my_app.solution[1])
	# heatmap(curl')
end

function main(;tstop=20.0, ntime=512, deltaRank=0, ml=1, cf=2, maxiter=10, fcf=1, relaxLyap=false, savegif=false, useTheta=false, k=4, ℜ=1600, deferDelta=(1, 1))
	my_app = KFlowApp(cf, useTheta, deltaRank, k, ℜ)

	tstart = 0.0

	core = XBraid.Init(
		comm, comm, tstart, tstop, ntime,
		my_step!, my_init, my_sum!, my_norm, my_access; app = my_app,
	)

	if deltaRank > 0
		XBraid.SetDeltaCorrection(core, deltaRank, my_basis_init, my_innerprod)
		# XBraid.SetDeferDelta(core, deferDelta...)
		XBraid.SetLyapunovEstimation(core; relax = relaxLyap, exponents = true)
	end
	XBraid.SetMaxLevels(core, ml)
	XBraid.SetMaxIter(core, maxiter)
	XBraid.SetCFactor(core, -1, cf)
	XBraid.SetNRelax(core, -1, fcf)
	XBraid.SetAccessLevel(core, 2)
	XBraid.SetAbsTol(core, 1e-6)

	XBraid.SetTimings(core, 2)
	XBraid.Drive(core)
	XBraid.PrintTimers(core)

	p = sortperm(my_app.times)

	if deltaRank > 0
	    exponents = sum(my_app.lyapunov_exps)
	    exponents = MPI.Allreduce(exponents, (+), comm)
	    exponents ./= tstop - tstart
	    if MPI.Comm_rank(comm) == 0
	        println("exponents: ", exponents)
	    end
	end

	# if deltaRank > 0
	# 	movingaverage(g, n) = [i < n ? mean(g[begin:i]) : mean(g[i-n+1:i]) for i in eachindex(g)]
	# 	println("Lyap Exps: ", sum(my_app.lyapunov_exps) ./ tstop)
	# 	exps = reduce(hcat, my_app.lyapunov_exps[p])'
	# 	nt = size(exps)[1]
	# 	Δt = tstop / nt
	# 	exps = movingaverage.([exps[:, i] ./ Δt for i ∈ 1:size(exps)[2]], length(exps[:, 1]))
	# 	exps_plot = plot(exps; legend = false)
	# 	savefig(exps_plot, "lyapunov_exps.png")
	# end

	if savegif
		preprocess(app, u) = spectralInterpolation(∇X(app, u), 128)'
		heatmapArgs = Dict(:ticks => false, :colorbar => false, :aspect_ratio => :equal)
		# anim = @animate for i in p[1:cf:end]
		count = 1
		for i in p
			u = my_app.solution[i]
			plots = [heatmap(preprocess(my_app, u); heatmapArgs...)]
			for j in 1:min(8, deltaRank)
				ψⱼ = my_app.lyapunov_vecs[i][j]
				push!(plots, heatmap(preprocess(my_app, ψⱼ); heatmapArgs...))
			end
			fig = plot(plots...; size = (600, 600))
			savefig(fig, "kflow_gif/beamer/kolmo_$(nₓ)_$(ntime)_ml$(ml)-$(count).png")
			count += 1
		end
		# gif(anim, "kflow_gif/kolmo_$(nₓ)_$(ntime)_ml$(ml).gif", fps = 20)
	end


	return my_app, core
end
