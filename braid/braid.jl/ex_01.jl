# main ex-01
include("XBraid.jl")
using .XBraid
using MPI

MPI.Init()
comm = MPI.COMM_WORLD

#=
using Ref{Float64} here because regular Float64
values are immutable and don't have stable addresses.
(Ref{Float64} is essentially a 1 element vector)
=#
function my_init(app, t)
    t == 0.0 && return Ref(1.0)
    return Ref(0.456)
end

# This must mutate u in place
function my_step!(app, status, u, ustop, tstart, tstop)
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
function my_access(app, status, u)
    XBraid.status_GetWrapperTest(status) && return
    t = XBraid.status_GetT(status)
    ti = XBraid.status_GetTIndex(status)
    print("t: $(t[]),\tu: $(u[])\n")
end

my_norm(app, u) = abs(u[])

test = false
if test
    test_app = XBraid.BraidApp(nothing, comm, my_step!, my_init, my_sum!, my_norm, my_access)
    # XBraid.postInitPrecompile(test_app)

    open("ex_01.test.out", "w") do file
        cfile = Libc.FILE(file)
        # XBraid.testInitAccess(test_app, 0.0, cfile)
        XBraid.testClone(test_app, 0.0, cfile)
        # XBraid.testSum(test_app, 0.0, cfile)
        # XBraid.testSpatialNorm(test_app, 0.0, cfile)
        XBraid.testBuf(test_app, 0.0, cfile)
    end
    MPI.Barrier(comm)
end

ntime = 10
tstart = 0.0
tstop = tstart + ntime / 2.0;
core = XBraid.Init(comm, comm, tstart, tstop, ntime, my_step!, my_init, my_sum!, my_norm, my_access)

XBraid.SetPrintLevel(core, 2)
XBraid.SetMaxLevels(core, 2)
XBraid.SetAbsTol(core, 1.e-6)
XBraid.SetCFactor(core, -1, 2)

XBraid.Warmup(core)
XBraid.SetTimings(core, 2)
XBraid.Drive(core; warmup=false)
XBraid.PrintTimers(core)
# no need for braid_Destroy(core), julia will take care of it :)
