using ForwardDiff, DiffResults, LinearAlgebra, PreallocationTools
using IterativeSolvers, LinearMaps, SparseArrays, Interpolations
using MPI, BenchmarkTools
using Plots

include("../../braid/braid.jl/XBraid.jl")
using .XBraid

MPI.Init()
comm = MPI.COMM_WORLD

Float = Float64

"""
struct to package preallocated caches for storing temp coords and velocity
as well as the sparse poisson matrix needed by project_incompressible!()
"""
struct TGApp
	x_d::DiffCache
	u_d::DiffCache
	ϕ_d::DiffCache
	Δ_s::SparseMatrixCSC
	solution
	lyapunov_vecs
	lyapunov_exps
	times
end

# default constructor
function TGApp(x::AbstractArray, u::AbstractArray, ϕ::AbstractArray, Δ::SparseMatrixCSC)
	TGApp(DiffCache(x), DiffCache(u), DiffCache(ϕ), Δ, [], [], [], [])
end

"""
reshapes 1D array into vector of 3D arrays, with no allocation
"""
function get_views(u)
	@views begin
		u_x = reshape(u[1:nₓ^3], (nₓ, nₓ, nₓ))
		u_y = reshape(u[nₓ^3+1:2*nₓ^3], (nₓ, nₓ, nₓ))
		u_z = reshape(u[2*nₓ^3+1:end], (nₓ, nₓ, nₓ))
	end
	return [u_x, u_y, u_z]
end

# some helpers for doing modulo arithmetic with CartesianIndex
Base.mod1(I::CartesianIndex, J::CartesianIndex) = CartesianIndex(broadcast(mod1, Tuple(I), Tuple(J)))
δ(I::CartesianIndex{N}, dim) where N = CartesianIndex(ntuple(i -> i == dim ? 1 : 0, N))

function TaylorGreen!(u, k = 1)
	kx = trunc(Int64, k/2) + k%2
	ky = trunc(Int64, k/2) + 1
	kz = trunc(Int64, k/2) + 1
	u_x, u_y, u_z = u
	x, y, z = coordinates
	@. u_x = sin(kx * x) * cos(ky * y) * cos(kz * z)
	@. u_y = -cos(kx * x) * sin(ky * y) * cos(kz * z)
	@. u_z = 0.0
	return
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
	@views @inbounds for I ∈ R
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
	out_arr = zeros(3 * nₓ^3)
	out = get_views(out_arr)
	Vx, Vy, Vz = V
	x, y, z = 1:3
	out[x] .= ∂(Vz, y) - ∂(Vy, z)
	out[y] .= ∂(Vx, z) - ∂(Vz, x)
	out[z] .= ∂(Vy, x) - ∂(Vx, y)
	return out
end

∇(F) = [∂(F, i) for i ∈ 1:3]
# Δ(F) = ∇_dot(∇(F)) # this is asymmetric...

function my_poisson(f_arr)
	out_arr = similar(f_arr)
	F = reshape(f_arr, (nₓ, nₓ, nₓ))
	out = reshape(out_arr, (nₓ, nₓ, nₓ))
	R = CartesianIndices(F)
	bound = last(R)
	@inbounds for I ∈ R
		out[I] = -6 * F[I]
		@inbounds for dim ∈ 1:3
			I⁺ = mod1(I + δ(I, dim), bound)
			I⁻ = mod1(I - δ(I, dim), bound)
			out[I] += F[I⁺] + F[I⁻]
		end
		out[I] *= 1 / Δx^2
	end
	# need a single point condition u(0,0,0) = 0., else
	# this system is singular (ones(nₓ^3) nullspace)
	out[first(R)] += F[first(R)]
	return out_arr
end

function my_poisson_mat(nₓ)
	flat_index(I::CartesianIndex) = I[1] + (I[2] - 1) * nₓ + (I[3] - 1) * nₓ^2
	rowinds = Array{Int}([])
	colinds = Array{Int}([])
	nzs = Array{Float}([])
	function make_connection!(I_r, I_c, nz)
		push!(rowinds, flat_index(I_r))
		push!(colinds, flat_index(I_c))
		push!(nzs, nz)
	end
	diag = -6 / Δx^2
	offd = 1 / Δx^2
	bound = CartesianIndex((nₓ, nₓ, nₓ))
	for I in CartesianIndices((1:nₓ, 1:nₓ, 1:nₓ))
		make_connection!(I, I, diag)
		for dim ∈ 1:3
			I⁺ = mod1(I + δ(I, dim), bound)
			I⁻ = mod1(I - δ(I, dim), bound)
			make_connection!(I, I⁺, offd)
			make_connection!(I, I⁻, offd)
		end
	end
	# point condition u(0,0,0) = 0., else system is singular
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
	x_range = range(0., 2π - Δx, nₓ)
	for i in 1:3
		itp = interpolate(u[i], BSpline(Linear(Periodic())))
		sitp = scale(itp, x_range, x_range, x_range)
		extp = extrapolate(sitp, Periodic())

		u_new[i] .= extp.(coords_new...)
	end
	u_arr .= u_new_arr
	return u_arr
end

"""
Solve Poisson problem for incompressible velocity
Δϕ = ∇⋅u
u = u - ∇ϕ
"""
function project_incompressible!(app::TGApp, u_arr)
	u = get_views(u_arr)
	ϕ_arr = get_tmp(app.ϕ_d, u_arr)
	ϕ_arr .= 0.
	rhs = vec(∇_dot(u))
	cg!(ϕ_arr, app.Δ_s, rhs; abstol=1/2*Δx^2)
	# cg!(ϕ_arr, app.Δ_s, rhs)

	ϕ = reshape(ϕ_arr, (nₓ, nₓ, nₓ))
	for i in 1:3
		u[i] .-= ∂(ϕ, i)
	end
	return u_arr
end

function base_step(app::TGApp, u_arr, Δt)
	# print("Advection : ")
	advect_semi_lagrangian!(app, u_arr, Δt)
	# diffuse_backward_Euler!(u_arr, Δt) # numerical diffusion may be enough
	# print("Projection: ")
	project_incompressible!(app, u_arr)
	return u_arr
end

nₓ = 64
Δx = 2π / nₓ

function my_init(app, t)
	u_arr = zeros(Float, 3 * nₓ^3)
	u = get_views(u_arr)
	TaylorGreen!(u)
	return u_arr
end

function my_basis_init(app, t, i)
	ψ_arr = zeros(Float, 3 * nₓ^3)
	ψ = get_views(ψ_arr)
	TaylorGreen!(ψ, i + 1)
	ψ_arr .+= Δx^3*randn(Float, 3 * nₓ^3)
	return ψ_arr
end

function my_step!(app, status, u, ustop, tstart, tstop)
	u = base_step(app, u, tstop - tstart)
	return
end

function my_step!(app, status, u, ustop, tstart, tstop, Ψ)
	Δt = tstop - tstart
	rank = length(Ψ)
	Ψ_new = reduce(hcat, Ψ)
	perturb(r) = base_step(app, u + r' * Ψ, Δt)

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

function my_access(app, status, u)
	index = XBraid.status_GetTIndex(status)
	lyap_exps = XBraid.status_GetLocalLyapExponents(status)
	lyap_vecs = XBraid.status_GetBasisVectors(status)
	if length(lyap_vecs) > 0
		push!(app.lyapunov_vecs, deepcopy(reduce(hcat, lyap_vecs)))
		push!(app.lyapunov_exps, lyap_exps)
	end

	push!(app.solution, deepcopy(u))
	push!(app.times, index)
end

my_norm(app, u) = Δx^3 * LinearAlgebra.norm2(u)
my_innerprod(app, u, v) = u' * v

coords_arr = zeros(Float, 3 * nₓ^3)
coordinates = get_views(coords_arr)
x_bound = 2π - Δx
x_range = range(0, x_bound, nₓ)
coordinates[1] .= [x for x ∈ x_range, y ∈ x_range, z ∈ x_range]
coordinates[2] .= [y for x ∈ x_range, y ∈ x_range, z ∈ x_range]
coordinates[3] .= [z for x ∈ x_range, y ∈ x_range, z ∈ x_range];

function test()
	# preallocations
	x_d = zeros(Float, 3 * nₓ^3)
	u_d = zeros(Float, 3 * nₓ^3)
	ϕ_d = zeros(Float, nₓ^3)
	Δ = my_poisson_mat(nₓ)
	my_app = TGApp(x_d, u_d, ϕ_d, Δ)

	# test user routines:
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

	curlz = ∇X(get_views(my_app.solution[1]))[3]
	heatmap(curlz[:, :, 1])
end

function main(;ml=1, tstop=4., ntime=32, delta_rank=0)
	# preallocations
	x_d = zeros(Float, 3 * nₓ^3)
	u_d = zeros(Float, 3 * nₓ^3)
	ϕ_d = zeros(Float, nₓ^3)
	Δ = my_poisson_mat(nₓ)
	my_app = TGApp(x_d, u_d, ϕ_d, Δ)

	# theme(:dracula)
	tstart = 0.0

	core = XBraid.Init(
		comm, comm, tstart, tstop, ntime,
		my_step!, my_init, my_sum!, my_norm, my_access; app = my_app,
	)::XBraid.BraidCore

	if delta_rank > 0
		XBraid.SetDeltaCorrection(core, delta_rank, my_basis_init, my_innerprod)
		XBraid.SetLyapunovEstimation(core; exponents=true)
	end
	XBraid.SetMaxLevels(core, ml)
	XBraid.SetMaxIter(core, 10)
	XBraid.SetCFactor(core, -1, 2)
	XBraid.SetAccessLevel(core, 1)
	XBraid.SetNRelax(core, -1, 0)
	XBraid.SetAbsTol(core, 1e-6)

	@time XBraid.Drive(core)

	p = sortperm(my_app.times)
	curlz_start = ∇X(get_views(my_app.solution[p][1]))[3]
	curlz_end = ∇X(get_views(my_app.solution[p][end]))[3]
	plot(heatmap(curlz_start[1, :, :]), heatmap(curlz_end[1, :, :]), size=(900,400))
	return my_app
end