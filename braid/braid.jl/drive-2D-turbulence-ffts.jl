using ForwardDiff, DiffResults, LinearAlgebra, PreallocationTools
using Interpolations, FFTW
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
struct KFlowApp
	cf::Int
	useTheta::Bool
    deltaRank::Int
    coords::Array{Float, 3}
	κ::Frequencies{Float}
	solution::Vector{Array{Float, 3}}
	lyapunov_vecs::Vector{Vector{Array{Float, 3}}}
	lyapunov_exps::Vector{Vector{Float}}
	times::Vector{Int}
	x_d::DiffCache{Array{Float, 3}}
	u_d::DiffCache{Array{Float, 3}}
	ϕ_d::DiffCache{Array{Complex{Float}, 2}}
    P̂::FFTW.FFTWPlan
end

# default constructor
function KFlowApp(cf::Integer, useTheta::Bool, deltaRank::Integer)
    coords = zeros(Float, nₓ, nₓ, 2)
    x_bound = 2π - Δx
    x_range = range(0, x_bound, nₓ)
    @views @inbounds for x ∈ x_range, y ∈ x_range 
        coords[x, y, 1] = x
        coords[x, y, 2] = y
    end
    κ = fftfreq(nₓ, nₓ)
	# preallocations
	x_d = zeros(Float, nₓ, nₓ, 2)
	u_d = zeros(Float, nₓ, nₓ, 2)
	ϕ_d = zeros(Float, nₓ, nₓ)
	KFlowApp(cf, useTheta, deltaRank, coords, κ, [], [], [], [], DiffCache(x_d), DiffCache(u_d), DiffCache(ϕ_d))
end

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
	ux, uy = u[:, :, 1], u[:, :, 2]
	kx, ky = get_wavenumber2D(i)
	fourierMode2D!(ux, kx)
	fourierMode2D!(uy, ky)
end

function TaylorGreen!(app::KFlowApp, u, k = 1)
    coords = app.coords
	kx = trunc(Int, k / 2) + k % 2
	ky = trunc(Int, k / 2) + 1
	@views u_x, u_y = u[:, :, 1], u[:, :, 3]
	@views x, y = coords[:, :, 1], coords[:, :, 2]
	@. u_x = sin(kx * x) * cos(ky * y)
	@. u_y = -cos(kx * x) * sin(ky * y)
	return
end

function kolmogorovForce!(u, Δt; k = 4, μ = 1.0)
	# Fₖ(x, y) = (0, sin(ky))
	@views @. u[:, :, 1] += μ * Δt * sin(k * coords[:, :, 2])
	return u
end

"""
compute self advection of u
"""
function advect_semi_lagrangian!(app::KFlowApp, u, Δt)
	# get cached arrays (either real or dual, depending on typeof(u))
	coords_new = get_tmp(app.x_d, u)
	u_new = get_tmp(app.u_d, u)

	# (1) trace each grid point backwards along u
	@. coords_new = coords - Δt * u

	# (2) interpolate the velocity field at the new grid points
	x_range = range(0.0, 2π - Δx, nₓ)
	for i in 1:2
		itp = interpolate(u[i], BSpline(Linear(Periodic())))
		sitp = scale(itp, x_range, x_range)
		extp = extrapolate(sitp, Periodic())

		u_new[i] .= @views extp.(coords_new[:, :, 1], coords_new[:, :, 2])
	end
	u .= u_new
	return u
end

"""
Solve Poisson problem for incompressible velocity field
Δϕ = ∇⋅u
u = u - ∇ϕ
"""
function project_incompressible!(app::KFlowApp, u::AbstractArray)
    κ = get_tmp(app.x_d, u)
	ϕ = get_tmp(app.ϕ_d, u)
    P̂ = app.P̂
    for κ_x ∈ app.κ, κ_y ∈ app.κ
        κ[x, y, 1] = κ_x
        κ[x, y, 2] = κ_y
    end
    û = P̂ * u
    @views κ_x, κ_y = κ[:, :, 1], κ[:, :, 2]
    @views û_x, û_y = û[:, :, 1], û[:, :, 2]
    # ϕ = inv(Δ)∇⋅u
    @. ϕ = -(1.0im*κ_x*û_x + 1.0im*κ_y*û_y)/(κ_x^2 + k_y^2)
    # u = u - ∇ϕ
    @. u -= [1.0im*κ_x*ϕ ;;; 1.0im*κ_y*ϕ]

	return u
end

function ∇X(app::KFlowApp, V)
    κ = get_tmp(app.x_d, u)
    for κ_x ∈ app.κ, κ_y ∈ app.κ
        κ[x, y, 1] = κ_x
        κ[x, y, 2] = κ_y
    end
	f = get_tmp(app.ϕ_d, u)
    V̂ = app.P̂ * V
	@views V̂x, V̂y = V̂[:, :, 1], V̂[:, :, 2]
	@views κx, κy = κ[:, :, 1], κ[:, :, 1]
	out = @. 1.0im*κx*V̂y - 1.0im*κy*V̂x
	return real(ifft(out))
end

function base_step!(app::KFlowApp, u::Vector{<: Number}, Δt::Float)
	kolmogorovForce!(u, Δt)
	advect_semi_lagrangian!(app, u, Δt)
	project_incompressible!(app, u)
	return u
end

function my_init(app::KFlowApp, t::Float)
	u = zeros(Float, nₓ, nₓ, 2)
	TaylorGreen!(app, u)
	return u
end


function my_basis_init(app::KFlowApp, t::Float, i::Integer)
	ψ = zeros(Float, nₓ, nₓ, 2)
	fourierMode2DVec!(ψ, i + 1)
	return ψ
end

function my_step!(
	app::KFlowApp, status::Ptr{Cvoid},
	u, ustop,
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
	app::KFlowApp, status::Ptr{Cvoid},
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

function my_access(app::KFlowApp, status, u)
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
			push!(app.lyapunov_vecs, deepcopy(vecs))
		end
	end
end

my_norm(app, u) = LinearAlgebra.normInf(u)
my_innerprod(app, u, v) = vec(u)' * vec(v)

const nₓ = 128
const Δx = 2π / nₓ

function test()

	my_app = KFlowApp(x_d, u_d, ϕ_d, Δ, 1, false)
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

	my_app = KFlowApp(x_d, u_d, ϕ_d, Δ, cf, useTheta)

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