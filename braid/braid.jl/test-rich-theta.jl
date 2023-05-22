include("XBraid.jl")
using .XBraid

using ToeplitzMatrices, LinearAlgebra
using Plots
using MPI
using IterativeSolvers

struct AdvDifApp
    nx::Int
    Δx::Float64
    A::Circulant
    I::SymmetricToeplitz
    useTheta::Bool
    useRich::Bool
    cf::Int
    order::Int
    solTf::Vector{Vector{Float64}}
    residuals::Vector{Float64}
end
# default constructor
function AdvDifApp(nx::Integer, Δx::Real, ν::Real, α::Real, useTheta::Bool, useRich::Bool, cf::Integer, order::Integer) 
    I = SymmetricToeplitz([1, zeros(nx-1)...])
    # Δ = (1/Δx^2)SymmetricToeplitz([-2, 1, zeros(nx-2)...])
    # ∇ = (1/Δx)Toeplitz([0, -1/2, zeros(nx-2)...], [0, 1/2, zeros(nx-2)...])
    Δ = (1/Δx^2)Circulant([-2, 1, zeros(nx-3)..., 1])
    ∇ = (1/Δx)Circulant([0, -1/2, zeros(nx-3)..., 1/2])
    # ∇ = (1/Δx)Circulant([-1, 1, zeros(nx-2)...])
    A = ν*Δ - α*∇
    AdvDifApp(nx, Δx, A, I, useTheta, useRich, cf, order, [], [])
end

function my_init(app, t)
    t == 0.0 && return sin.(range(0., 2π-app.Δx, app.nx))
    return randn(app.nx)
end

function backward_euler!(app, u, Δt)
    u .= (app.I - Δt * app.A) \ u
    # u .= cg(app.I - Δt * app.A, u)
    return u
end

function θ_euler!(app, u, Δt, m)
    θ = (m - 1)/2m
    u .= (app.I - (1-θ)Δt * app.A) \ ((app.I + θ*Δt * app.A) * u)
    return u
end

function SDIRK2!(app, u, Δt)
    a = 1/√2
    k1 = Δt*((app.I - (1-a) * Δt * app.A) \ (app.A*u))
    k2 = Δt*((app.I - (1-a) * Δt * app.A) \ (app.A*(u + (2a-1)*k1)))
    u .+= 1/2*(k1 + k2)
    return u
end

function θSDIRK2!(app, u, Δt, m)
    α = √2√(m*(m-1))/2m
    θ = (m^2 -3m - √2√(m*(m-1)))/(2m^2 - 4m)
    k1 = Δt*((app.I - (1-α) * Δt * app.A) \ (app.A*u))
    k2 = Δt*((app.I - (1-α) * Δt * app.A) \ (app.A*(u + (2α-1)*k1)))
    u .+= (θ)k1 + (1-θ)k2
    return u
end

function base_step!(app, u, Δt, level)
    if app.order == 2
        if app.useTheta
            throw("Theta method not implemented for order 2")
        end
        return SDIRK2!(app, u, Δt)
    end

    if app.useTheta && level > 0
        return θSDIRK2!(app, u, Δt, app.cf^level)
        # return θ_euler!(app, u, Δt, app.cf^level)
    end

    return backward_euler!(app, u, Δt)
end

function my_step!(app, status, u, ustop, tstart, tstop)
    Δt = tstop - tstart
    level = XBraid.status_GetLevel(status)
    if !app.useRich || level == 0
        return base_step!(app, u, Δt, level)
    end
    m = app.cf ^ level
    p = app.order
    if app.useTheta
        p += 1
    end
    θ = 2^p * (m^p - 1) / (m^p * (2^p - 1))
    u_sub = deepcopy(u)
    base_step!(app, u_sub, Δt/2, level)
    base_step!(app, u_sub, Δt/2, level)
    base_step!(app, u, Δt, level)
    @. u = θ*u_sub + (1-θ)*u
    return u
end

function my_access(app, status, u)
    XBraid.status_GetWrapperTest(status) && return
    ti = XBraid.status_GetTIndex(status)
    ntime = XBraid.status_GetNTPoints(status)
    if ti == ntime
        push!(app.solTf, deepcopy(u))
    end
end

function test(;ν=1., α=1., order=1)
    MPI.Init()
    comm = MPI.COMM_WORLD

    my_app = AdvDifApp(512, 2π/513, ν, α, false, false, 2, order)
    app = XBraid.BraidApp(my_app, comm, my_step!, my_init, my_access)

    XBraid.testInitAccess(app, 0.)
    XBraid.testSpatialNorm(app, 0.)
    u = my_init(my_app, 0.)
    plot(u, label="u0", legend=:bottomright)
    my_step!(my_app, nothing, u, nothing, 0., .1)
    plot!(u, label="u1")
end


function main(;tstop=2π, ntime=512, nx=512, ν=1., α=1., useTheta=false, useRich=false, cf=2, ml=2, order=1)
    MPI.Init()
    comm = MPI.COMM_WORLD

    Δx = 2π / nx
    app = AdvDifApp(nx, Δx, ν, α, useTheta, useRich, cf, order)

    core = XBraid.Init(comm, 0.0, tstop, ntime, my_step!, my_init, my_access; app=app)

    XBraid.SetPrintLevel(core, 2)
    XBraid.SetMaxLevels(core, ml)
    XBraid.SetCFactor(core, -1, app.cf)
    XBraid.SetAccessLevel(core, 1)
    XBraid.SetNRelax(core, -1, 1)
    XBraid.SetAbsTol(core, 1e-8)
    XBraid.SetMaxIter(core, 40)
    XBraid.SetSkip(core, false)
    
    XBraid.Drive(core)

    append!(app.residuals, XBraid.GetRNorms(core))

    return app
end