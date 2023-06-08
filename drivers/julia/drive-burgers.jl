using LinearAlgebra, PreallocationTools, ForwardDiff, DiffResults
using MPI, Interpolations, IterativeSolvers
using SparseArrays, ToeplitzMatrices
using FFTW
using Plots
using LinearAlgebra: norm2
using Random: seed!
# theme(:dracula)

include("../../braid/braid.jl/XBraid.jl")
using .XBraid

MPI.Init()
comm = MPI.COMM_WORLD

struct BurgerApp
    x::Vector{Float64}
    κ::Frequencies{Float64}
    μ::Float64
    cf::Int
    useTheta::Bool
    # preallocated caches that support ForwardDiff:
    x_d::DiffCache
    y_d::DiffCache
    # storage for solution values
    solution::Vector{Vector{Float64}}
    lyap_vecs::Vector{Union{Matrix{Float64}, Vector{Float64}}}
    lyap_exps::Vector{Vector{Float64}}
    times::Vector{Int}
    # pre-planned fourier transform
    P̂
end

function BurgerApp(x::Vector{<:Real}, κ::Frequencies{<:Real}, μ::Real, cf::Integer, useTheta::Bool, x_cache::AbstractArray, y_cache::AbstractArray)
    BurgerApp(x, κ, μ, cf, useTheta, DiffCache(x_cache), DiffCache(y_cache), [], [], [], [], plan_fft(x))
end

# system parameters (globally scoped)
const lengthScale = 2π
const nₓ = 1024
const Δx = lengthScale / nₓ

function stencil_to_circulant(stencil::Vector{T}, nx::Integer) where T
    @assert isodd(length(stencil)) "stencil should have odd number of entries"
    stencilWidth = length(stencil)
    center = ceil(Int, length(stencil)/2)
    columnView = @view(stencil[end:-1:1])
    v = [columnView[center:end]; zeros(T, nx-stencilWidth); columnView[1:center-1]]
    return Circulant(v)
end

function euler_mat(Δt::Float64; μ::Float64=.00001)
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

function semi_lagrangian!(burger::BurgerApp, y, v, Δt)
    x_back = get_tmp(burger.x_d, y)
    y_intp = get_tmp(burger.y_d, y)

    # initialize interpolation
    itp = interpolate(y, BSpline(Linear(Periodic())))
    sitp = scale(itp, range(0., lengthScale - Δx, nₓ))
    extp = extrapolate(sitp, Periodic())

    x_back .= burger.x .- Δt * v
    
    # y_intp .= interp_periodic.(x_back, Ref(y))
    y_intp .= extp.(x_back)
    y .= y_intp
    return y
end

function diffuse_beuler!(burger::BurgerApp, y, Δt::Float64; init_guess=nothing)
    y_tmp = get_tmp(burger.y_d, y)
    y_tmp .= y

    if init_guess !== nothing
        y .= init_guess
    end

    cg!(y, euler_mat(Δt), y_tmp)
    return y
end

function diffuse_fft!(burger::BurgerApp, y::AbstractArray, Δt::Real)
    P̂ = burger.P̂
    κ = burger.κ
    # y .= real(P̂ \ (exp.(Δt*(κ.^2 - κ.^4))))
    ŷ = P̂ * y
    @. ŷ *= exp(-burger.μ * Δt * κ^2) # standard diffusion
    # @. ŷ *= exp(Δt*(κ^2 - κ^4)) # KS-equation operator
    y .= real(P̂ \ ŷ)
    return y
end

function extractPartials(y::Vector{ForwardDiff.Dual{T,V,P}}) where {T,V,P}
    ps = zeros(V, size(y)..., P)
    for i ∈ eachindex(y), j ∈ 1:P
        @inbounds ps[i, j] = ForwardDiff.partials(y[i], j)
    end
    return ps
end

function fillDualArray!(y::Vector{ForwardDiff.Dual{T,V,P}}, vs, ps) where {T,V,P}
    checkbounds(ps, firstindex(y), 1:P)
    for i ∈ eachindex(y)
        @inbounds y[i] = ForwardDiff.Dual{T}(vs[i], ntuple(j -> @inbounds(ps[i, j]), P))
    end
end

# this enables ForwardDiff through the FFT where it normally doesn't work
function diffuse_fft!(burger::BurgerApp, y::Vector{ForwardDiff.Dual{T,V,P}}, Δt::Real) where {T, V, P}
    vs = ForwardDiff.value.(y)
    ps = extractPartials(y)
    diffuse_fft!(burger, vs, Δt)
    map(eachcol(ps)) do p diffuse_fft!(burger, p, Δt) end
    fillDualArray!(y, vs, ps)
    return y
end

# user routines:
function my_init(burger::BurgerApp, t::Float64)
    # sin
    u = sin.(2π/lengthScale * burger.x)
    # gaussian
    # u = exp.(-(burger.x .- lengthScale/2).^2 / (lengthScale/5)^2)
    # if t == 0
    #     seed!(1234)
    #     u .+= 1e-1randn(nₓ) .* u
    # end
    # @. u = exp(-(burger.x - lengthScale/2)^2)
    return u
end

function my_basis_init(burger::BurgerApp, t::Float64, k::Int32)
    x = burger.x
    ψ = similar(x)
    if k % 2 == 0
        @. ψ = cos(k/2*x)
    else
        @. ψ = sin((k+1)/2*x)
    end
    return ψ
end

function base_step!(burger::BurgerApp, u, Δt::Float64; init_guess=nothing)
    u_mid = deepcopy(u)
    # semi_lagrangian!(burger, u_mid, u, Δt/2)
    # semi_lagrangian!(burger, u, u_mid, Δt) # u(∇ u)
    # semi_lagrangian!(burger, u, sin.(burger.x), Δt) # ∇ u
    semi_lagrangian!(burger, u, u, Δt) # u(∇ u)
    if burger.μ > 0.
        diffuse_fft!(burger, u, Δt) # Δu
    end
end

function my_step!(
    burger::BurgerApp, status::XBraid.Status.StepStatus,
    u, ustop, tstart::Float64, tstop::Float64
)
    Δt = tstop - tstart
    level = XBraid.Status.getLevel(status)
    if !burger.useTheta || level == 0
        base_step!(burger, u, Δt; init_guess=ustop)
    else
        # richardson based θ method
        m = burger.cf^level
        θ = 2(m - 1)/m
        # p = 0.86
        # p = 1.
        # θ = 2^p * (m^p - 1) / (m^p * (2^p - 1))
        # θ = 2
        u_sub = deepcopy(u)
        base_step!(burger, u_sub, Δt/2; init_guess=(ustop .- u)./2)
        base_step!(burger, u_sub, Δt/2; init_guess=ustop)
        base_step!(burger, u, Δt; init_guess=u_sub)
        @. u = θ*u_sub + (1-θ)*u
    end
    return u
end

function my_step!(burger::BurgerApp, status::XBraid.Status.StepStatus, u, ustop, tstart::Float64, tstop::Float64, Ψ)
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

function my_access(burger::BurgerApp, status::XBraid.Status.AccessStatus, u)
    XBraid.Status.getWrapperTest(status) && return
    ti = XBraid.Status.getTIndex(status)
    push!(burger.solution, deepcopy(u))
    push!(burger.times, ti)
    if XBraid.Status.getDeltaRank(status) > 0
        Ψ = XBraid.Status.getBasisVectors(status)
        λ = XBraid.Status.getLocalLyapExponents(status)
        push!(burger.lyap_vecs, deepcopy(reduce(hcat, Ψ)))
        push!(burger.lyap_exps, deepcopy(λ))
    end
end

# test user routines:
function test()
    x_new = zeros(nₓ)
    y_new = zeros(nₓ)
    x = Array(range(0., lengthScale-Δx, nₓ))
    κ = 2π/lengthScale .* fftfreq(nₓ, nₓ)
    burger = BurgerApp(x, κ, 1e-3, 4, false, x_new, y_new);
    test_app = XBraid.BraidApp(burger, comm, my_step!, my_init, my_access; basis_init=my_basis_init)

    open("drive-burgers.test.out", "w") do file
        cfile = Libc.FILE(file)
        XBraid.testBuf(test_app, 0.0, cfile)
        XBraid.testSpatialNorm(test_app, 0.0, cfile)
        XBraid.testDelta(test_app, 0.0, 0.1, 3, cfile)
    end
    # plot!(burger.lyap_vecs[1][1])
    # plot(burger.solution[1])
end

function main(;tstop=π, ntime=nₓ, deltaRank=0, useTheta=false, fmg=false, ml=1, cf=4, saveGif=false, maxiter=30, μ=0.)
    x_new = zeros(nₓ)
    y_new = zeros(nₓ)
    x = Array(range(0., lengthScale-Δx, nₓ))
    κ = 2π/lengthScale .* fftfreq(nₓ, nₓ)
    burger = BurgerApp(x, κ, μ, 4, useTheta, x_new, y_new);

    tstart = 0.0
    core = XBraid.Init(comm, tstart, tstop, ntime, my_step!, my_init, my_access; app=burger)

    if deltaRank > 0
        XBraid.setDeltaCorrection(core, deltaRank, my_basis_init)
        XBraid.setLyapunovEstimation(core; exponents=true)
    end

    # println("Wrapper test:")
    # @time test()

    XBraid.setMaxIter(core, maxiter)
    XBraid.setMaxLevels(core, ml)
    XBraid.setCFactor(core, -1, cf)
    XBraid.setAccessLevel(core, 1)
    XBraid.setNRelax(core, -1, 1)
    XBraid.setAbsTol(core, 1e-6)
    XBraid.setSkip(core, false)
    if fmg
        XBraid.setFMG(core)
        XBraid.setNFMG(core, 1)
    end

    XBraid.setTimings(core, 2)
    XBraid.Drive(core)
    XBraid.printTimers(core)

    if deltaRank > 0
        exponents = sum(burger.lyap_exps)
        exponents = MPI.Allreduce(exponents, (+), comm)
        exponents ./= tstop
        if MPI.Comm_rank(comm) == 0
            println("exponents: ", exponents)
        end
    end

    if saveGif
        for (ti, u) in zip(burger.times, burger.solution)
            plt = plot(burger.x, u; ylim=(-1.5, 1.5))
            savefig(plt, "burgers_gif/$(lpad(ti, 6, "0")).png")
        end
        MPI.Barrier(comm)

        if MPI.Comm_rank(comm) == 0
            println("animating...")
            fnames = [lpad(i, 6, "0") for i in 0:ntime]
            anim = Animation("burgers_gif", fnames)
            Plots.buildanimation(anim, "burgers_gif/anim.gif", fps=30)
        end
    end

    return burger, XBraid.getRNorms(core)
end

# test()
# main();
nothing