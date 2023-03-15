using LinearAlgebra, PreallocationTools, ForwardDiff, DiffResults
using MPI, Interpolations, IterativeSolvers
using SparseArrays, ToeplitzMatrices
using BenchmarkTools, Plots
theme(:dark)

include("XBraid.jl")
using .XBraid

MPI.Init()
comm = MPI.COMM_WORLD


struct burger_app
    x::Vector{Float64}
    cf::Integer
    # preallocated caches that support ForwardDiff:
    x_d::DiffCache
    y_d::DiffCache
    solution
    lyap_vecs
    lyap_exps
    times
end
function burger_app(x::Vector{Float64}, cf::Integer, x_d::DiffCache, y_d::DiffCache)
    burger_app(x, cf, x_d, y_d, [], [], [], [])
end

# system parameters (globally scoped)
lengthScale = 2π
nₓ = 128
Δx = lengthScale / nₓ
η = .2
x = Array(range(0.0, lengthScale - Δx, nₓ))

function stencil_to_circulant(stencil::Vector{T}, nx::Integer) where T
    @assert isodd(length(stencil)) "stencil should have odd number of entries"
    stencilWidth = length(stencil)
    center = ceil(Int, length(stencil)/2)
    columnView = reverse(stencil)
    v = [columnView[center:end]; zeros(T, nx-stencilWidth); columnView[1:center-1]]
    return Circulant(v)
end

function euler_mat(Δt; μ=.00001)
    # diffusion
    # A = spdiagm(-1 => ones(nₓ-1), 0 => -2*ones(nₓ), 1 => ones(nₓ-1))
    # A[1, end] = 1
    # A[end, 1] = 1
    # A = sparse(I, (nₓ, nₓ)) - Δt/Δx^2 * μ * A

    # KS operator
    ∇⁴ = [1, -4, 6, -4, 1]
    ∇² = [0, 1, -2, 1, 0]
    I  = [0, 0, 1, 0, 0]
    stencil = @. I - Δt*(-1/Δx^4 * ∇⁴ - 1/Δx^2 * ∇²)
    # stencil = @. I - Δt(-μ/Δx^4 * ∇⁴)
    stencil_to_circulant(stencil, nₓ)
end

function semi_lagrangian!(burger::burger_app, y, Δt)
    x_back = get_tmp(burger.x_d, y)
    y_intp = get_tmp(burger.y_d, y)

    x_back .= x .- Δt * y
    itp = interpolate(y, BSpline(Linear(Periodic())))
    sitp = scale(itp, range(0., lengthScale - Δx, nₓ))
    extp = extrapolate(sitp, Periodic())
    # y_intp .= interp_periodic.(x_back, Ref(y))
    y_intp .= extp.(x_back)
    y .= y_intp
    return y
end

function diffuse_beuler!(burger::burger_app, y, Δt; init_guess=nothing)
    y_tmp = get_tmp(burger.y_d, y)
    y_tmp .= y

    if init_guess !== nothing
        y .= init_guess
    end

    cg!(y, euler_mat(Δt), y_tmp)
    return y
end

# preallocations (performance critical, passed as arguments)
x_new = zeros(nₓ)
y_new = zeros(nₓ)

# user routines:
function my_init(burger, t)
    u = similar(burger.x)
    # u .= sin.(2π/lengthScale * burger.x)
    @. u = exp(-(x - lengthScale/2)^2)
    return u
end

function my_basis_init(burger, t, k)
    x = burger.x
    ψ = similar(x)
    if k % 2 == 0
        @. ψ = cos(k/2*x)
    else
        @. ψ = sin((k+1)/2*x)
    end
    return ψ
end

function base_step!(burger, u, Δt; init_guess=nothing)
    # diffuse_beuler!(burger, u, Δt/2; init_guess=init_guess)
    u = semi_lagrangian!(burger, u, Δt)
    # diffuse_beuler!(burger, u, Δt; init_guess=init_guess)
end

function my_step!(burger, status, u, ustop, tstart, tstop)
    Δt = tstop - tstart
    level = XBraid.status_GetLevel(status)
    # if level == 0
    if true
        base_step!(burger, u, Δt; init_guess=ustop)
    else
        # richardson based θ method
        m = burger.cf^level
        θ = 2(m - 1)/m
        # θ = 2
        u_sub = deepcopy(u)
        base_step!(burger, u_sub, Δt/2; init_guess=(ustop .- u)./2)
        base_step!(burger, u_sub, Δt/2; init_guess=ustop)
        base_step!(burger, u, Δt; init_guess=u_sub)
        @. u = θ*u_sub + (1-θ)*u
    end
    return u
end

function my_step!(burger, status, u, ustop, tstart, tstop, Ψ)
    rank = length(Ψ)
    Ψ_new = reduce(hcat, Ψ)
    # perturb(r) = base_step!(burger, u + r' * Ψ, Δt)
    perturb(r) = my_step!(burger, status, u + r' * Ψ, ustop, tstart, tstop)

    result = DiffResults.DiffResult(u, Ψ_new)
    result = ForwardDiff.jacobian!(result, perturb, zeros(rank))
    for i in eachindex(Ψ)
        Ψ[i] .= Ψ_new[:, i]
    end
end

function my_sum!(burger, a, x, b, y)
    @. y = a*x + b*y
end

function my_access(burger, status, u)
    push!(burger.solution, deepcopy(u))
    t = XBraid.status_GetT(status)
    ti = XBraid.status_GetTIndex(status)
    push!(burger.times, t)
    Ψ = XBraid.status_GetBasisVectors(status)
    push!(burger.lyap_vecs, deepcopy(Ψ))
    λ = XBraid.status_GetLocalLyapExponents(status)
    push!(burger.lyap_exps, deepcopy(λ))
end

my_norm(burger, u) = LinearAlgebra.norm2(u)
my_innerprod(burger, u, v) = u' * v


# test user routines:
function test()
    burger = burger_app(x, 4, DiffCache(x_new), DiffCache(y_new));
    test_app = XBraid.BraidApp(
        burger, comm, comm,
        my_step!, my_init,
        my_sum!, my_norm, my_access,
        my_basis_init, my_innerprod)

    XBraid.testBuf(test_app, 0.0)
    XBraid.testSpatialNorm(test_app, 0.0)
    XBraid.testDelta(test_app, 0.0, 0.1, 3)
    plot(burger.solution[1])
    plot!(burger.lyap_vecs[1][1])
end

function main(;tstop=5., ntime=256, delta_rank=0, ml=1, cf=4)
    tstart = 0.0
    burger = burger_app(x, cf, DiffCache(x_new), DiffCache(y_new));

    core = XBraid.Init(
        comm, comm, tstart, tstop, ntime,
        my_step!, my_init, my_sum!, my_norm, my_access; app=burger
    )

    if delta_rank > 0
        XBraid.SetDeltaCorrection(core, delta_rank, my_basis_init, my_innerprod)
        XBraid.SetLyapunovEstimation(core; exponents=true)
    end
    XBraid.SetMaxIter(core, 45)
    XBraid.SetMaxLevels(core, ml)
    XBraid.SetCFactor(core, -1, cf)
    XBraid.SetAccessLevel(core, 1)
    XBraid.SetNRelax(core, -1, 1)
    XBraid.SetAbsTol(core, 1e-6)

    @time XBraid.Drive(core)

    if delta_rank > 0
        print(sum(burger.lyap_exps)/tstop)
    end

    p = sortperm(burger.times)
    plot(burger.solution[p][1:64:end])

    return burger
end