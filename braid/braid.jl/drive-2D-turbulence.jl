using ForwardDiff, DiffResults, LinearAlgebra, PreallocationTools
using IterativeSolvers, SparseArrays, Interpolations
using MPI, BenchmarkTools, Plots, LaTeXStrings
using Statistics: mean

include("XBraid.jl")
using .XBraid

MPI.Init()
comm = MPI.COMM_WORLD

Float = Float64

"""
struct to package preallocated caches for storing temp coords and velocity
as well as the sparse poisson matrix needed by project_incompressible!
"""
struct TGApp
	x_d::DiffCache{Vector{Float}}
	u_d::DiffCache{Vector{Float}}
	ϕ_d::DiffCache{Vector{Float}}
	Δ_s::SparseMatrixCSC{Float, Int}
	cf::Int
	useTheta::Bool
	solution::Vector{Vector{Float}}
	lyapunov_vecs::Vector{Matrix{Float}}
	lyapunov_exps::Vector{Vector{Float}}
	times::Vector{Int}
end

# default constructor
function TGApp(x::AbstractArray, u::AbstractArray, ϕ::AbstractArray, Δ::SparseMatrixCSC, cf::Integer, useTheta::Bool)
	TGApp(DiffCache(x), DiffCache(u), DiffCache(ϕ), Δ, cf, useTheta, [], [], [], [])
end

"""
reshapes 1D array into vector of 2D arrays, with no allocation
"""
function get_views(u)
	@views begin
		u_x = reshape(u[1:nₓ^2], (nₓ, nₓ))
		u_y = reshape(u[nₓ^2+1:2*nₓ^2], (nₓ, nₓ))
	end
	return [u_x, u_y]
end

# some helpers for doing modulo arithmetic with CartesianIndex
Base.mod1(I::CartesianIndex, J::CartesianIndex) = CartesianIndex(broadcast(mod1, Tuple(I), Tuple(J)))
δ(I::CartesianIndex{N}, dim) where N = CartesianIndex(ntuple(i -> i == dim ? 1 : 0, N))

oneball(n) = [(n + 1 - i, i) for i ∈ 1:n]
wavenumbers(shell) = reduce(hcat, [oneball(i) for i ∈ 1:shell])
is_in_shell(i, shell) = i <= trunc(Int, shell * (shell + 1) / 2)
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

function fourierMode1D(x::Real, k::Integer)
	k % 2 == 1 && return cos(trunc(k / 2) * x)
	k % 2 == 0 && return sin(trunc((k + 1) / 2) * x)
end

function fourierMode2D!(u_field, i)
	kx, ky = get_wavenumber2D(i)
	x, y = coords
	@. u_field = fourierMode1D(x, kx) * fourierMode1D(y, ky)
end

function fourierMode2DVec!(u, i)
	ux, uy = u
	kx, ky = get_wavenumber2D(i)
	fourierMode2D!(ux, kx)
	fourierMode2D!(uy, ky)
end

function TaylorGreen!(u, k = 1)
	kx = trunc(Int, k / 2) + k % 2
	ky = trunc(Int, k / 2) + 1
	u_x, u_y = u
	x, y = coords
	@. u_x = sin(kx * x) * cos(ky * y)
	@. u_y = -cos(kx * x) * sin(ky * y)
	return
end

function kolmogorovForce!(u, Δt; k = 4, μ = 1.0)
	# Fₖ(x, y) = (0, sin(ky))
	@. @views u[1] += μ * Δt * sin(k * coords[2])
	return u
end

"""
finite difference approx first derivative
of A wrt dimension *dim*
"""
function ∂(A, dim; h = Δx)
	# This works for A of arbitrary dimension
	out = similar(A)
	R = CartesianIndices(A)
	bound = last(R)
	for I ∈ R
		# second order
		I⁻ = mod1(I - δ(I, dim), bound)
		I⁺ = mod1(I + δ(I, dim), bound)
		out[I] = 1 / 2h * (A[I⁺] - A[I⁻])
	end
	return out
end

# divergence, curl, gradient, and laplacian
∇_dot(V) =
	sum(enumerate(V)) do (i, V_i)
		∂(V_i, i)
	end
function ∇X(V)
	out = zeros(nₓ, nₓ)
	Vx, Vy = V
	x, y = 1:2
	out .= ∂(Vy, x) - ∂(Vx, y)
	return out
end
∇(F) = [∂(F, i) for i ∈ 1:3]
# Δ(F) = ∇_dot(∇(F)) # this is asymmetric...

function my_poisson_mat(nₓ)
	flat_index(I::CartesianIndex) = I[1] + (I[2] - 1) * nₓ
	rowinds = Array{Int}([])
	colinds = Array{Int}([])
	nzs = Array{Float}([])
	function make_connection!(I_r, I_c, nz)
		push!(rowinds, flat_index(I_r))
		push!(colinds, flat_index(I_c))
		push!(nzs, nz)
	end
	diag = -4 / Δx^2
	offd = 1 / Δx^2
	bound = CartesianIndex((nₓ, nₓ))
	for I in CartesianIndices((1:nₓ, 1:nₓ))
		make_connection!(I, I, diag)
		for dim ∈ 1:2
			I⁺ = mod1(I + δ(I, dim), bound)
			I⁻ = mod1(I - δ(I, dim), bound)
			make_connection!(I, I⁺, offd)
			make_connection!(I, I⁻, offd)
		end
	end
	# point condition u(0,0) = 0., else system is singular
	nzs[1] += 1.0
	return sparse(rowinds, colinds, nzs)
end

"""
compute self advection of u
"""
function advect_semi_lagrangian!(app::TGApp, u_arr, Δt)
	# get cached arrays (either real or dual, depending on typeof(u))
	coords_new_arr = get_tmp(app.x_d, u_arr)
	u_new_arr = get_tmp(app.u_d, u_arr)

	u = get_views(u_arr)
	u_new = get_views(u_new_arr)

	# (1) trace each grid point backwards along u
	@. coords_new_arr = coords_arr - Δt * u_arr
	coords_new = get_views(coords_new_arr)

	# (2) interpolate the velocity field at the new grid points
	x_range = range(0.0, 2π - Δx, nₓ)
	for i in 1:2
		itp = interpolate(u[i], BSpline(Linear(Periodic())))
		sitp = scale(itp, x_range, x_range)
		extp = extrapolate(sitp, Periodic())

		u_new[i] .= extp.(coords_new...)
	end
	u_arr .= u_new_arr
	return u_arr
end

"""
Solve Poisson problem for incompressible velocity field
Δϕ = ∇⋅u
u = u - ∇ϕ
"""
function project_incompressible!(app::TGApp, u_arr::Vector{<: Number})
	u = get_views(u_arr)
	ϕ_arr = get_tmp(app.ϕ_d, u_arr)
	ϕ_arr .= 0.0
	rhs = reshape(∇_dot(u), nₓ^2)
	cg!(ϕ_arr, app.Δ_s, rhs; reltol = 1 / 2 * Δx^2)
	# cg!(ϕ_arr, app.Δ_s, rhs)
	ϕ = reshape(ϕ_arr, (nₓ, nₓ))
	for i in 1:2
		u[i] .-= ∂(ϕ, i)
	end
	return u_arr
end

function base_step!(app::TGApp, u_arr::Vector{<: Number}, Δt::Float)
	kolmogorovForce!(get_views(u_arr), Δt)
	advect_semi_lagrangian!(app, u_arr, Δt)
	# diffuse_backward_Euler!(u_arr, Δt) # numerical diffusion may be enough
	project_incompressible!(app, u_arr)
	return u_arr
end


function my_init(app::TGApp, t::Float)
	u_arr = zeros(Float, 2 * nₓ^2)
	u = get_views(u_arr)
	TaylorGreen!(u; k=3)
	# kolmogorovForce!(u, 1.; μ=1.)
	# u_arr .+= 1e-2*randn(Float, 2*nₓ^2)
	return u_arr
end


function my_basis_init(app::TGApp, t::Float, i::Integer)
	ψ_arr = zeros(Float, 2 * nₓ^2)
	# ψ_arr = randn(Float, 2*nₓ^2)
	ψ = get_views(ψ_arr)
	fourierMode2DVec!(ψ, i + 1)
	# ψ_arr += 1e-1*randn(Float, 2*nₓ^2)
	# TaylorGreen!(ψ, i + 1)
	return ψ_arr
end

function my_step!(
	app::TGApp, status::Ptr{Cvoid},
	u::Vector{<: Number}, ustop::Vector{<: Number},
	tstart::Float, tstop::Float,
)
	Δt = tstop - tstart
	level = XBraid.status_GetLevel(status)
	if !app.useTheta || level == 0
		u = base_step!(app, u, Δt)
	else
		# richardson based θ method
		m = app.cf^level
		θ = 2(m - 1) / m
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
	app::TGApp, status::Ptr{Cvoid},
	u::Vector{Float}, ustop::Vector{Float},
	tstart::Float, tstop::Float,
	Ψ::Vector{Any},
)
	rank = length(Ψ)
	Ψ_new = reduce(hcat, Ψ)
	perturb(r) = my_step!(app, status, u + r' * Ψ, ustop, tstart, tstop)

	result = DiffResults.DiffResult(u, Ψ_new)
	result = ForwardDiff.jacobian!(result, perturb, zeros(rank))
	for i in 1:rank
		Ψ[i] .= Ψ_new[:, i]
	end
	return
end

function my_sum!(app, a, x, b, y)
	@. y = a * x + b * y
end

function my_access(app::TGApp, status, u)
	XBraid.status_GetWrapperTest(status) && return
	index = XBraid.status_GetTIndex(status)
	if index % app.cf == 0
		push!(app.solution, deepcopy(u))
		push!(app.times, index)
		rank = XBraid.status_GetDeltaRank(status)
		if rank > 0
			exps = XBraid.status_GetLocalLyapExponents(status)
			vecs = XBraid.status_GetBasisVectors(status)
			push!(app.lyapunov_exps, exps)
			push!(app.lyapunov_vecs, deepcopy(reduce(hcat, vecs)))
		end
	end
end

my_norm(app, u) = LinearAlgebra.normInf(u)
my_innerprod(app, u, v) = u' * v

nₓ = 128
Δx = 2π / nₓ

coords_arr = zeros(Float, 2 * nₓ^2)
coords = get_views(coords_arr)
x_bound = 2π - Δx
x_range = range(0, x_bound, nₓ)
coords[1] .= [x for x ∈ x_range, y ∈ x_range]
coords[2] .= [y for x ∈ x_range, y ∈ x_range]

function test()
	# preallocations
	x_d = zeros(Float, 2 * nₓ^2)
	u_d = zeros(Float, 2 * nₓ^2)
	ϕ_d = zeros(Float, nₓ^2)
	Δ = my_poisson_mat(nₓ)

	my_app = TGApp(x_d, u_d, ϕ_d, Δ, 1, false)
	test_app = XBraid.BraidApp(
		my_app, comm, comm,
		my_step!, my_init,
		my_sum!, my_norm, my_access,
		my_basis_init, my_innerprod)

	XBraid.testInitAccess(test_app, 0.0)
	XBraid.testClone(test_app, 0.0)
	XBraid.testSpatialNorm(test_app, 0.0)
	XBraid.testBuf(test_app, 0.0)
	@time XBraid.testDelta(test_app, 0.0, 0.1, 3)

	curl = ∇X(get_views(my_app.solution[1]))
	heatmap(curl')
end

function main(; tstop = 64.0, ntime = 512, deltaRank = 0, ml = 1, cf = 2, maxiter = 10, fcf = 1, relaxLyap = false, savegif = true, useTheta = false, deferDelta = (1, 1))
	# preallocations
	x_d = zeros(Float, 2 * nₓ^2)
	u_d = zeros(Float, 2 * nₓ^2)
	ϕ_d = zeros(Float, nₓ^2)
	Δ = my_poisson_mat(nₓ)

	my_app = TGApp(x_d, u_d, ϕ_d, Δ, cf, useTheta)

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
	XBraid.SetAccessLevel(core, 1)
	XBraid.SetAbsTol(core, 1e-6)

	XBraid.SetTimings(core, 2)
	XBraid.Drive(core)
	XBraid.PrintTimers(core)

	p = sortperm(my_app.times)

	if savegif
		heatmapArgs = Dict(:ticks => false, :colorbar => false, :aspect_ratio => :equal, :c => :plasma)
		anim = @animate for i in p
			u = get_views(my_app.solution[i])
			# plots = [heatmap(∇_dot(u)'; heatmapArgs...)]
			plots = [heatmap(∇X(u)'; heatmapArgs...)]
			for j in 1:min(8, deltaRank)
				ψⱼ = get_views(my_app.lyapunov_vecs[i][:, j])
				# push!(plots, heatmap(∇_dot(ψⱼ)'; heatmapArgs...))
				push!(plots, heatmap(∇X(ψⱼ)'; heatmapArgs...))
			end
			plot(plots...; size = (600, 600))
		end
		gif(anim, "kflow_gif/kolmo_$(nₓ)_$(ntime)_ml$(ml).gif", fps = 20)
	end

    # if deltaRank > 0
    #     exponents = sum(my_app.lyapunov_exps)
    #     exponents = MPI.Allreduce(exponents, (+), comm)
    #     exponents ./= tstop
    #     if MPI.Comm_rank(comm) == 0
    #         println("exponents: ", exponents)
    #     end
    # end

	if deltaRank > 0
		movingaverage(g, n) = [i < n ? mean(g[begin:i]) : mean(g[i-n+1:i]) for i in eachindex(g)]
		println("Lyap Exps: ", sum(my_app.lyapunov_exps) ./ tstop)
		exps = reduce(hcat, my_app.lyapunov_exps[p])'
		nt = size(exps)[1]
		Δt = tstop / nt
		exps = movingaverage.([exps[:, i] ./ Δt for i ∈ 1:size(exps)[2]], length(exps[:, 1]))
		exps_plot = plot(exps; legend = false)
		savefig(exps_plot, "lyapunov_exps.png")
	end

	return my_app
end
