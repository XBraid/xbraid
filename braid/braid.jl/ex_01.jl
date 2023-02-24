# main ex-01
include("XBraid.jl")
using .XBraid
using ForwardDiff, MPI

MPI.Init()
comm = MPI.COMM_WORLD

#=
using Ref(Float64) here because regular Float64
values are immutible and don't have stable addresses.
e.g.

julia> function f(x)
            x = 2
        end
f (generic function with 1 method)

julia> x = 4
4

julia> f(x)
2

julia> x
4
=#
function my_init(app, t)
    t == 0.0 && return Ref(1.0)
    return Ref(0.456)
end

# This must mutate u in place
function my_step!(app, u, ustop, tstart, tstop)
    # backward Euler
    u[] = 1.0 / (1.0 + tstop - tstart) * u[]
end

#=
    compute the sum
    y ← ax + by
    this must mutate y in place
=#
function my_sum!(app, a, x, b, y)
    y[] = a * x[] + b * y[]
end

#=
    called by the solver to give access to 
    solution values.
=#
function my_access(app, u)
    println(u[])
end

my_norm(app, u) = abs(u[])

test_app = XBraid.BraidApp(nothing, comm, my_step, my_init, my_sum, my_norm, my_access)

XBraid.testInitAccess(test_app, 0.0)
XBraid.testClone(test_app, 0.0)
XBraid.testSum(test_app, 0.0)
XBraid.testSpatialNorm(test_app, 0.0)
XBraid.testBuf(test_app, 0.0)

ntime = 10
tstart = 0.0
tstop = tstart + ntime / 2.0;
core = XBraid.Init(comm, comm, tstart, tstop, ntime, my_step, my_init, my_sum, my_norm, my_access)

XBraid.SetPrintLevel(core, 2)
XBraid.SetMaxLevels(core, 2)
XBraid.SetAbsTol(core, 1.e-6)
XBraid.SetCFactor(core, -1, 2)

XBraid.Drive(core)
# no need for braid_Destroy(core), julia will take care of it :)
