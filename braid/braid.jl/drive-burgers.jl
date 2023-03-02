using LinearAlgebra, PreallocationTools, ForwardDiff, DiffResults
using MPI

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
const nₓ = 64
const Δx = 2π / nₓ
const Δt = 0.1
x = Array(range(0.0, 2π - Δx, nₓ))

# preallocations (performance critical, passed as arguments)
x_new = zeros(nₓ)
y_new = zeros(nₓ)

# user routines:
function interp_periodic(coord, f)
    bound = nₓ
    scaled = coord / Δx
    i = trunc(scaled)
    r = scaled - i
    i = mod1(Int(i + 1), bound)
    (1 - r) * f[i] + r * f[mod1(i + 1, bound)]
end

function semi_lagrangian!(burger::burger_app, y, Δt)
    x_back = get_tmp(burger.x_d, y)
    y_intp = get_tmp(burger.y_d, y)

    x_back .= x .- Δt * y
    y_intp .= interp_periodic.(x_back, Ref(y))
    y .= y_intp
    return y
end

function my_init(burger, t)
    u = similar(burger.x)
    u .= sin.(x) .+ 1e-2*randn(length(x))
    # u .= sin.(x)
    return u
end

function my_basis_init(burger, t, k)
    ψ = similar(burger.x)
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
    return
end

function my_step!(burger, status, u, ustop, tstart, tstop, Ψ)
    Δt = tstop - tstart
    rank = length(Ψ)
    Ψ_new = reduce(hcat, Ψ)
    # perturb(r) = semi_lagrangian!(burger, u + Ψ_new * r, Δt)
    perturb(r) = semi_lagrangian!(burger, u + r' * Ψ, Δt)

    result = DiffResults.DiffResult(u, Ψ_new)
    result = ForwardDiff.jacobian!(result, perturb, zeros(rank))
    for i in 1:length(Ψ)
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
    XBraid.testDelta(test_app, 0.0, Δt, delta_rank)
end

burger = burger_app(x, DiffCache(x_new), DiffCache(y_new));
tstart = 0.0
tstop = 4.0
ntime = 128
delta_rank = 1

core = XBraid.Init(
    comm, comm, tstart, tstop, ntime,
    my_step!, my_init, my_sum!, my_norm, my_access; app=burger
)

XBraid.SetDeltaCorrection(core, delta_rank, my_basis_init, my_innerprod)
XBraid.SetLyapunovEstimation(core, false, true, true)
XBraid.SetMaxLevels(core, 3)
XBraid.SetCFactor(core, -1, 2)
XBraid.SetAccessLevel(core, 1)
XBraid.SetNRelax(core, -1, 0)
XBraid.Drive(core)
