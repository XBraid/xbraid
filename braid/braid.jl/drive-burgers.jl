using LinearAlgebra, PreallocationTools, ForwardDiff, DiffResults
using MPI, SparseArrays, Interpolations, IterativeSolvers
using BenchmarkTools, Plots
unicodeplots()
theme(:dark)

include("XBraid.jl")
using .XBraid

MPI.Init()
comm = MPI.COMM_WORLD


struct burger_app
    x
    # preallocated caches that support ForwardDiff:
    x_d::DiffCache
    y_d::DiffCache
    solution # stores the total solution vector over space-time
    lyap_vecs # stores the lyapunov vectors
    lyap_exps # stores the lyapunov exponents
    times
end
burger_app(x, x_d, y_d) = burger_app(x, x_d, y_d, [], [], [], [])

# system parameters (globally scoped)
nₓ = 128
Δx = 2π / nₓ
Δt = 0.1
η = .2
x = Array(range(0.0, 2π - Δx, nₓ))

euler_mat = spdiagm(-1 => ones(nₓ-1), 0 => -2*ones(nₓ), 1 => ones(nₓ-1))
euler_mat[1, end] = 1
euler_mat[end, 1] = 1
euler_mat = spdiagm(ones(nₓ)) - Δt * η * euler_mat

# preallocations (performance critical, passed as arguments)
x_new = zeros(nₓ)
y_new = zeros(nₓ)

# user routines:
function interp_periodic(coord, f)
    f_circ = CircularArray(f)
    scaled = coord / Δx
    i = trunc(Int64, scaled)
    r = scaled - i
    i = i + 1
    (1 - r) * f_circ[i] + r * f_circ[i + 1]
end

function semi_lagrangian!(burger::burger_app, y, Δt)
    x_back = get_tmp(burger.x_d, y)
    y_intp = get_tmp(burger.y_d, y)

    x_back .= x .- Δt * y
    itp = interpolate(y, BSpline(Linear(Periodic())))
    sitp = scale(itp, range(0., 2π - Δx, nₓ))
    extp = extrapolate(sitp, Periodic())
    # y_intp .= interp_periodic.(x_back, Ref(y))
    y_intp .= extp.(x_back)
    y .= y_intp
    return y
end

function diffuse_beuler!(burger::burger_app, y, Δt)
    y_tmp = get_tmp(burger.y_d, y)
    y_tmp .= y
    cg!(y, euler_mat, y_tmp)
    return y
end

function my_init(burger, t)
    u = similar(burger.x)
    u .= sin.(burger.x)
    # u .= sin.(x)
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

function my_step!(burger, status, u, ustop, tstart, tstop)
    Δt = tstop - tstart
    u = semi_lagrangian!(burger, u, Δt)
    u = diffuse_beuler!(burger, u, Δt)
    return
end

function my_step!(burger, status, u, ustop, tstart, tstop, Ψ)
    Δt = tstop - tstart
    rank = length(Ψ)
    Ψ_new = reduce(hcat, Ψ)
    perturb(r) = diffuse_beuler!(burger, semi_lagrangian!(burger, u + r' * Ψ, Δt), Δt)

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
if true
    burger = burger_app(x, DiffCache(x_new), DiffCache(y_new));
    test_app = XBraid.BraidApp(
        burger, comm, comm,
        my_step!, my_init,
        my_sum!, my_norm, my_access,
        my_basis_init, my_innerprod)

    XBraid.testBuf(test_app, 0.0)
    XBraid.testSpatialNorm(test_app, 0.0)
    XBraid.testDelta(test_app, 0.0, Δt, 3)
end

burger = burger_app(x, DiffCache(x_new), DiffCache(y_new));
tstart = 0.0
tstop = 4.0
ntime = 256
delta_rank = 1

core = XBraid.Init(
    comm, comm, tstart, tstop, ntime,
    my_step!, my_init, my_sum!, my_norm, my_access; app=burger
)

# XBraid.SetDeltaCorrection(core, delta_rank, my_basis_init, my_innerprod)
# XBraid.SetLyapunovEstimation(core, false, true, true)
XBraid.SetMaxLevels(core, 3)
XBraid.SetCFactor(core, -1, 4)
XBraid.SetAccessLevel(core, 1)
XBraid.SetNRelax(core, -1, 1)
XBraid.SetAbsTol(core, 1e-6)

@time XBraid.Drive(core)

p = sortperm(burger.times)
plot(burger.solution[p][1:64:end])
