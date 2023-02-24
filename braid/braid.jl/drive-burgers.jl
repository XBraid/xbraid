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
end

function interp_periodic(coord, f)
    bound = nₓ
    # r, i = modf(coord / Δx)
    scaled = coord / Δx
    i = trunc(Int, scaled)
    r = scaled - i
    i = mod1(i + 1, bound)
    (1 - r) * f[i] + r * f[mod1(i + 1, bound)]
end

function semi_lagrangian(burger::burger_app, y, Δt)
    # load preallocated arrays from caches 
    # (loads normal data for y::Real, diffcache for y::Dual)
    x_back = get_tmp(burger.x_d, y)
    y_intp = get_tmp(burger.y_d, y)

    x_back .= x .- Δt * y
    y_intp .= interp_periodic.(x_back, Ref(y))
    y .= y_intp
    return y
end

# provide function which allocates (important) a new vector
# and initializes it.
function my_init(burger, t)
    # similar(x) allocates a new vector with the shape of x,
    # but with uninitialized entries
    u = similar(x)
    @. u = sin(x)
    return u
end

# initialize a basis vector at time t and spatial index k ∈ [0, 1, ..., nₓ-1]
# here using the real fourier modes
function my_basis_init(burger, t, k)
    ψ = similar(burger.x)
    if k % 2 == 0
        @. ψ = cos(k / 2 * x)
    else
        @. ψ = sin((k + 1) / 2 * x)
    end
    return ψ
end

# provide two functions, one to take a normal time-step...
function my_step!(burger, u, ustop, tstart, tstop)
    Δt = tstop - tstart
    u = semi_lagrangian(burger, u, Δt)
    return
end

# ... and one to take a time-step and propagate Ψ
function my_step!(burger, u, ustop, tstart, tstop, Ψ)
    rank = length(Ψ)
    Δt = tstop - tstart
    # Ψ will be a vector of arrays... not ideal, but that's how XBraid stores them
    # so we just concatenate the vectors together into a 2D array
    Ψ_new = reduce(hcat, Ψ)
    # differentiating this function wrt the Ψ coordinate vector r
    # yields the propagated lyapunov vectors (∂ᵤΦ) * Ψ
    perturb(r) = semi_lagrangian(burger, u + r' * Ψ, Δt)

    result = DiffResults.DiffResult(u, Ψ_new)
    result = ForwardDiff.jacobian!(result, perturb, zeros(rank))
    for i in 1:rank
        Ψ[i] .= Ψ_new[:, i]
    end
end

# functions for sum, norm, and inner product of vectors
function my_sum!(burger, a, x, b, y)
    # make sure we overwrite y in place
    @. y = a * x + b * y
end

function my_access(burger, u)
    # println(u[nₓ÷4])
end

my_norm(burger, u) = LinearAlgebra.norm2(u)
my_innerprod(burger, u, v) = u' * v

# system parameters (globally scoped)
const nₓ = 64
const Δx = 2π / nₓ
const Δt = 0.1
x = Array(range(0.0, 2π - Δx, nₓ))

# preallocations (performance critical, passed as arguments)
x_new = zeros(nₓ)
y_new = zeros(nₓ)

burger = burger_app(
    x,
    DiffCache(x_new),
    DiffCache(y_new)
);

# test user routines:
if false
    test_app = XBraid.BraidApp(
        burger, comm, comm,
        my_step!, my_init,
        my_sum!, my_norm, my_access,
        my_basis_init, my_innerprod)

    XBraid.testBuf(test_app, 0.0)
    XBraid.testSpatialNorm(test_app, 0.0)
    XBraid.testDelta(test_app, 0.0, Δt, delta_rank)
end

tstart = 0.0
tstop = 5.0
ntime = 512
delta_rank = 32

core = XBraid.Init(
    comm, comm, tstart, tstop, ntime,
    my_step!, my_init, my_sum!, my_norm, my_access; app=burger
)
# XBraid.SetDeltaCorrection(core, delta_rank, my_basis_init, my_innerprod)
XBraid.SetMaxLevels(core, 2)
XBraid.SetCFactor(core, -1, 4)
XBraid.SetAccessLevel(core, 1)

XBraid.Drive(core)
