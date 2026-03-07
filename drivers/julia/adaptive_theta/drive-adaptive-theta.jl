using ToeplitzMatrices, LinearAlgebra, IterativeSolvers
using Plots
using MPI

include("../../../braid/braid.jl/XBraid.jl")
include("butcher-tables.jl")
using .XBraid

mutable struct AdvDifApp
    # usual parameters for simulation
    nx::Int
    Δx::Float64
    A::Circulant
    I::SymmetricToeplitz
    useTheta::Bool
    useRich::Bool
    order::Int
    cf::Int
    max_levels::Int
    solTf::Vector{Vector{Float64}}
    residuals::Vector{Float64}
    times::Vector{Float64}
    spatial_tol::Float64

    # parameters for adaptive theta
    skip_downcycle::Bool
    cgorder::Int
    flag_refine_downcycle::Bool
    fine_btable::ButcherTable
    cg_btables::Vector{Vector{ButcherTable}} # butcher tables for each level/step
end

fine_btables = [beuler, sdirk212, esdirk3, sdirk4]

# default constructor
function AdvDifApp(nx::Integer, Δx::Real, ν::Real, α::Real, useTheta::Bool, useRich::Bool, order::Integer, cf::Integer, ml::Integer; spatial_tol=1e-3, cgorder=order+1, skip=true)
    I = SymmetricToeplitz([1, zeros(nx-1)...])
    Δ = (1/Δx^2)Circulant([-2, 1, zeros(nx-3)..., 1])
    ∇ = (1/Δx)Circulant([0, -1/2, zeros(nx-3)..., 1/2])
    # ∇ = (1/Δx)Circulant([-1, 1, zeros(nx-2)...])
    A = ν*Δ - α*∇
    AdvDifApp(nx, Δx, A, I, useTheta, useRich, order, cf, ml, [], [], [], spatial_tol, skip, cgorder, false, fine_btables[order], [])
end

mutable struct ThetaVector
    u::Vector{Float64}
    ψ::Vector{Float64}
    t_prior::Float64
end
function ThetaVector(u, order::Integer)
    numconditions = [0, 1, 3, 7, 16]
    order > 5 && return ThetaVector(u, zeros(numconditions[5]), 0.)
    return ThetaVector(u, zeros(numconditions[order]), 0.)
end

function dirk_step!(app::AdvDifApp, u::Vector{<:Real}, Δt::Real, btable::ButcherTable; embedding=nothing)
    k = zeros(app.nx, btable.s)
    for i in 1:btable.s
        k[:, i] .= Δt * app.A * (u .+ sum(btable.A[i, j]*k[:, j] for j ∈ 1:i-1; init=zeros(app.nx)))
        if btable.A[i, i] != 0
            k[:, i] .= (app.I -  Δt * btable.A[i, i] * app.A) \ k[:, i]
        end
    end
    u .= u + sum(btable.b[j] * k[:, j] for j ∈ 1:btable.s)
    if embedding !== nothing
        @assert length(embedding) == btable.s
        err_est = sum(embedding[j] * k[:, j] for j ∈ 1:btable.s) - sum(btable.b[j] * k[:, j] for j ∈ 1:btable.s)
        return norm(err_est)
    end
    return 0.
end

function get_btable(app::AdvDifApp, status::XBraid.Status.StepStatus, ti::Integer)
    level = XBraid.Status.getLevel(status)
    if level == 0 || !app.useTheta
        return app.fine_btable
    end
    iu, il = XBraid.Status.getTIUL(status, level)
    i = ti - il + 1
    return app.cg_btables[level][i]
end

#=====================
XBraid functions
=====================#

function my_init(app, t)
    t == 0.0 && return ThetaVector(sin.(range(0., 2π-app.Δx, app.nx)), app.cgorder)
    u = randn(app.nx)
    # u = zeros(app.nx)
    return ThetaVector(u, app.cgorder)
end

T = forest
function my_step!(app::AdvDifApp, status::XBraid.Status.StepStatus, u::ThetaVector, ustop::ThetaVector, tstart, tstop)
    Δt = tstop - tstart
    ti = XBraid.Status.getTIndex(status)
    caller = XBraid.Status.getCallingFunction(status)
    level = XBraid.Status.getLevel(status)

    # get Butcher table for this step
    btable = get_btable(app, status, ti)
    
    # f-relax
    if app.flag_refine_downcycle && caller ∈ [XBraid.CallerFRestrict, XBraid.CallerFASResidual]
        cfactor = XBraid.Status.getCFactor(status)
        # first step of f-interval
        if XBraid.isCPoint(ti, cfactor)
            u.ψ .= 0.
            u.t_prior = tstart
        end
        Δt_prior = tstart - u.t_prior
        ψ_old = deepcopy(u.ψ)
        ψ_step = [ψ(btable, t) for t ∈ T.keys[1:length(u.ψ)]]

        # second order
        u.ψ[T[:{τ}]] += Δt^2 * ψ_step[T[:{τ}]] + Δt * Δt_prior

        # third order
        if app.cgorder >= 3
            u.ψ[T[:{{τ}}]] += Δt^3 * ψ_step[T[:{{τ}}]] + Δt * ψ_old[T[:{τ}]] + Δt^2 * Δt_prior * ψ_step[T[:{τ}]]
            u.ψ[T[:{τ, τ}]] += Δt^3 * ψ_step[T[:{τ, τ}]] + 2Δt^2 * Δt_prior * ψ_step[T[:{τ}]] + Δt * Δt_prior^2
        end

        # last step of f-interval
        if XBraid.isCPoint(ti+1, cfactor)
            # normalize order conditions
            Δt_c = tstop - u.t_prior
            # 2nd order
            u.ψ[T[:{τ}]] /= Δt_c^2
            f!, J!, fill, guess = θsdirk2_lhs!, θsdirk2_J!, θsdirk2_fill, θsdirk2_guess
            # 3rd order
            if app.cgorder >= 3
                u.ψ[T[:{{τ}}]] /= Δt_c^3
                u.ψ[T[:{τ, τ}]] /= Δt_c^3
                f!, J!, fill, guess = θsdirk3_lhs!, θsdirk3_J!, θsdirk3_fill, θsdirk3_guess
            end

            # solve order conditions:
            # @info """
            # Solving order conditions at t = $tstop
            # ψ = $(u.ψ)
            # """
            iu, il = XBraid.Status.getTIUL(status, level+1)
            ic = XBraid.mapFineToCoarse(ti, cfactor) - il + 1
            app.cg_btables[level+1][ic] = solve_order_conditions(f!, J!, fill, app.cgorder, u.ψ; guess=guess)
        end
    end

    # turn off order conditions once upcycle starts
    if app.flag_refine_downcycle && caller == XBraid.Status.CallerFInterp
        @info "done computing coarse grid butcher tables"
        app.flag_refine_downcycle = false
    end

    # finally, take the step
    embed = nothing
    if level == 0 && app.order == 2
        embed = sdirk212_emb
    end
    # println("level: $level, ti: $ti")
    # display(btable)
    # println("emb: $embed")
    err_est = dirk_step!(app, u.u, Δt, btable; embedding=embed)
    if err_est >= app.spatial_tol && caller != XBraid.CallerWarmup
        XBraid.Status.setRFactor(status, 2)
    end
end

function my_sync!(app::AdvDifApp, status::XBraid.SyncStatus)
    app.useTheta || return
    caller = XBraid.Status.getCallingFunction(status)
    iter = XBraid.Status.getIter(status)
    # println("sync! caller: $caller, iter: $(XBraid.Status.getIter(status))")

    if caller == XBraid.CallerDrive_AfterInit || caller == XBraid.Status.CallerFRefine_AfterInitHier
        app.cg_btables = Vector{ButcherTable}[]
        # initialize storage for coarse grid butcher tables
        nlevels = XBraid.Status.getNLevels(status)
        for l ∈ 1:nlevels-1
            iu, il = XBraid.Status.getTIUL(status, l)
            ncpoints = iu - il + 1
            push!(app.cg_btables, [app.fine_btable for i ∈ 1:ncpoints])
        end

        # turn on computation of coarse grid butcher tables
        app.flag_refine_downcycle = true
    end

    # make sure we still compute these when we skip the first downcycle
    if app.skip_downcycle && caller == XBraid.CallerDrive_TopCycle && iter == 0
        app.flag_refine_downcycle = true
    end
    app.flag_refine_downcycle && @info "Recomputing coarse grid butcher tables"
    return
end

function my_norm(app::AdvDifApp, u::ThetaVector)
    LinearAlgebra.normInf(u.u)
end

function my_sum!(app::AdvDifApp, α::Real, x::ThetaVector, β::Real, y::ThetaVector)
    y.u .= α * x.u + β * y.u
end

function my_access(app::AdvDifApp, status, u)
    XBraid.Status.getWrapperTest(status) && return
    t = XBraid.Status.getT(status)
    ti = XBraid.Status.getTIndex(status)
    ntime = XBraid.Status.getNTPoints(status)
    done = XBraid.Status.getDone(status)
    if ti == ntime
        push!(app.solTf, deepcopy(u.u))
    end
    if done
        push!(app.times, t)
    end
end

function test(;ν=1., α=1., order=1)
    MPI.Init()
    comm = MPI.COMM_WORLD

    my_app = AdvDifApp(512, 2π/513, ν, α, false, false, order, 2, 2)
    app = XBraid.BraidApp(my_app, comm, my_step!, my_init, my_access; sum=my_sum!, spatialnorm=my_norm)

    XBraid.testInitAccess(app, 0.)
    XBraid.testSum(app, 0.)
    XBraid.testSpatialNorm(app, 0.)
    u = my_init(my_app, 0.)
    println(my_norm(my_app, u))
    plot(u.u, label="u0", legend=:bottomright)
    my_step!(my_app, XBraid.StepStatus(), u, u, 0., .05)
    println(my_norm(my_app, u))
    plot!(u.u, label="u1")
end

function main(;tstop=2π, ntime=512, nx=512, ν=1., α=1., useTheta=false, useRich=false, cf=4, ml=2, maxiter=10, order=1, skip=false, refine=false, maxrefine=2, finefcf=true, tol=1e-3)
    MPI.Init()
    comm = MPI.COMM_WORLD

    Δx = 2π / nx
    app = AdvDifApp(nx, Δx, ν, α, useTheta, useRich, order, cf, ml; skip=skip, spatial_tol=tol)

    core = XBraid.Init(comm, 0., tstop, ntime, my_step!, my_init, my_access;
                       app=app, sync=my_sync!, sum=my_sum!, spatialnorm=my_norm)

    XBraid.setPrintLevel(core, 2)
    XBraid.setMaxLevels(core, ml)
    XBraid.setCFactor(core, -1, app.cf)
    XBraid.setAccessLevel(core, 1)
    XBraid.setNRelax(core, -1, 1)
    XBraid.setAbsTol(core, 1e-8)
    XBraid.setMaxIter(core, maxiter)
    XBraid.setSkip(core, skip)
    XBraid.setRefine(core, refine)
    XBraid.setMaxRefinements(core, maxrefine)
    if !finefcf
        XBraid.setNRelax(core, 0, 0)
    end
    
    XBraid.Drive(core)

    append!(app.residuals, XBraid.getRNorms(core))

    return app
end