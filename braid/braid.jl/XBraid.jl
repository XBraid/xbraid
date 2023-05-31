#=BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory.
 * 
 * This file is part of XBraid. For support, post issues to the XBraid Github page.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER=#

module XBraid
# Not sure we want to export anything...
# export Init, Drive, etc.

using Serialization: serialize, deserialize, Serializer # used to pack arbitrary julia objects into buffers
using LinearAlgebra: norm2, dot
using MPI

include("BraidUtils.jl")
using .BraidUtils

#=
 |  For this to work, XBraid must be compiled as a shared library.
 |  $ make braid shared=yes
 |  and the Fortran flags braid_Fortran_SpatialCoarsen, braid_Fortran_Residual,
 |  braid_Fortran_TimeGrid, and braid_Fortran_Sync must all be set to zero in
 |  braid.h before compiling. (TODO: figure out why this is...)
 =#

# This can contain anything (TODO: do I even need this?)
mutable struct BraidVector{T}
    user_vector::T
    VecType::Type
end
BraidVector(u::T) where {T} = BraidVector(u, T)

OptionalFunction = Union{Function,Nothing}

"""
This is an internal structure that is not generally exposed to the user.
This enables the user interface to only use memory safe function calls.
User defined data structures that are not time-dependent can be declared in 
the global scope, or they can be packed into a user defined object and passed to
Init(), in which case the object will be passed as the first argument in step, sum, spatialnorm, and access.
_app.user_app is passed as the first argument to some user defined function.
Stores any time-independent data the user may need to compute a time-step.
Large, preallocate arrays should be packed into this object, since operations
involving globally scoped variables are slower than local variables.
struct braid_App end

This may be manually contructed in order to pass to test functions.

julia> app = BraidApp(my_app, MPI.COMM_WORLD, MPI.COMM_WORLD, my_step, my_init, my_sum, my_norm, my_access);

julia> testSpatialNorm(app, 0.);
 """
mutable struct BraidApp
    user_app::Any     # user defined app data structure

    comm::MPI.Comm    # global mpi communicator
    comm_t::MPI.Comm  # temporal mpi communicator

    # user functions
    step::Function
    init::Function
    sum::Function
    spatialnorm::Function
    access::Function

    sync::OptionalFunction
    basis_init::OptionalFunction
    inner_prod::OptionalFunction

    # dictionary to store globally scoped references to allocated braid_vector objects
    ref_ids::IdDict{UInt64,BraidVector}
    bufsize::Integer  # expected serialized size of user's vector
    bufsize_lyap::Integer
    user_AppType::Type
    user_VecType::Type
    user_BasType::Type
end

# default constructor
function BraidApp(app, comm::MPI.Comm, comm_t::MPI.Comm, step::Function, init::Function, sum::Function, norm::Function, access::Function, sync::OptionalFunction, basis_init::OptionalFunction, inner_prod::OptionalFunction)
    BraidApp(app, comm, comm_t, step, init, sum, norm, access, sync, basis_init, inner_prod, IdDict(), 0, 0, Nothing, Nothing, Nothing)
end

function BraidApp(app, comm::MPI.Comm, step, init, access; comm_t=comm, sum=default_sum!, spatialnorm=default_norm, sync=nothing, basis_init=nothing, inner_prod=default_inner_prod)
    BraidApp(app, comm, comm_t, step, init, sum, spatialnorm, access, sync, basis_init, inner_prod)
end



"""
Stores all the information needed to run XBraid. Create this object with Init(), then pass this to Drive() to run XBraid.

julia> core = XBraid.Init(MPI.COMM_WORLD, MPI.COMM_WORLD, 

julia> XBraid.Drive(core);

"""
mutable struct BraidCore
    # internal values
    _braid_core::Ptr{Cvoid}
    _braid_app::BraidApp
    tstart::Float64
    tstop::Float64
    ntime::Int32
    function BraidCore(_braid_core, _braid_app, tstart, tstop, ntime)
        x = new(_braid_core, _braid_app, tstart, tstop, ntime)
        finalizer(x) do core
            # @async println("Destroying BraidCore $(core._braid_core)")
            @ccall libbraid.braid_Destroy(core._braid_core::Ptr{Cvoid})::Cint
        end
    end
end


# import wrapper functions and status structures
include("Status.jl")
include("Wrappers.jl")
using .Status, .Wrappers

function postInitPrecompile(app::BraidApp)
    #= 
    	# This is a hack to make sure every processor knows how big the user's vector will be.
    	# We can also take this time to precompile the user's step function so it doesn't happen
    	# in serial on the coarse grid.
    	=#
    GC.enable(false) # disable garbage collection
    app_ptr = pointer_from_objref(app)::Ptr{Cvoid}
    pp = malloc_null_double_ptr(Cvoid)

    _jl_init!(app_ptr, 0.0, pp)
    u_ptr = unsafe_load(pp)
    u = unsafe_pointer_to_objref(u_ptr)::BraidVector
    VecType = typeof(u.user_vector)
    AppType = typeof(app.user_app)

    !precompile(app.step, (AppType, Status.StepStatus, VecType, VecType, Float64, Float64)) && println("failed to precompile step")
    !precompile(app.spatialnorm, (AppType, VecType)) && println("failed to precompile norm")
    !precompile(app.sum, (AppType, Float64, VecType, Float64, VecType)) && println("failed to precompile sum")
    if app.access !== nothing
        !precompile(app.access, (AppType, Status.AccessStatus, VecType)) && println("failed to precompile access")
    end

    # some julia functions that can be precompiled
    precompile(deepcopy, (BraidVector{VecType},))
    precompile(unsafe_store!, (Ptr{Float64}, Float64,))
    precompile(Tuple{typeof(Base.getproperty),BraidVector{VecType},Symbol})
    precompile(Tuple{typeof(Base.unsafe_store!),Ptr{Int32},Int64})
    precompile(Tuple{typeof(deserialize),Serializer{Base.GenericIOBuffer{Array{UInt8,1}}},DataType})
    precompile(Tuple{typeof(Base.unsafe_wrap),Type{Array{UInt8,1}},Ptr{UInt8},Int64})
    precompile(Tuple{Type{NamedTuple{(:read, :write, :maxsize),T} where T<:Tuple},Tuple{Bool,Bool,Int64}})

    _jl_free!(app_ptr, u_ptr)
    Base.Libc.free(pp)
    GC.enable(true) # re-enable garbage collection

    app.user_AppType = AppType
    app.user_VecType = VecType
end

function deltaPrecompile(app::BraidApp)
    GC.enable(false)
    app_ptr = pointer_from_objref(app)::Ptr{Cvoid}
    pp = malloc_null_double_ptr(Cvoid)
    AppType = app.user_AppType
    VecType = app.user_VecType

    _jl_init_basis!(app_ptr, 0.0, Int32(0), pp)
    ψ_ptr = unsafe_load(pp)
    ψ = unsafe_pointer_to_objref(ψ_ptr)
    BasType = typeof(ψ.user_vector)
    !precompile(app.step, (AppType, Status.StepStatus, VecType, VecType, Float64, Float64, Vector{BasType})) && println("failed to compile step_du")
    !precompile(app.inner_prod, (AppType, VecType, VecType)) && println("failed to compile inner_prod u⋅u")
    !precompile(app.inner_prod, (AppType, BasType, VecType)) && println("failed to compile inner_prod ψ⋅u")
    !precompile(app.inner_prod, (AppType, VecType, BasType)) && println("failed to compile inner_prod u⋅ψ")
    !precompile(app.inner_prod, (AppType, BasType, BasType)) && println("failed to compile inner_prod ψ⋅ψ")

    _jl_free!(app_ptr, ψ_ptr)
    Base.Libc.free(pp)
    GC.enable(true)

    app.user_BasType = BasType
end

"""
Create a new BraidCore object. The BraidCore object is used to run the XBraid solver, and is destroyed when the object is garbage collected.
"""
function Init(
    comm_world::MPI.Comm,
    tstart::Real,
    tstop::Real,
    ntime::Integer,
    step::Function,
    init::Function,
    access::Function;
    sum              = default_sum!,
    spatialnorm      = default_norm,
    sync             = nothing,
    comm_t::MPI.Comm = comm_world,
    app              = nothing
)::BraidCore
    _app = BraidApp(app, comm_world, comm_t, step, init, sum, spatialnorm, access, sync, nothing, nothing)
    _core_ptr = malloc_null_double_ptr(Cvoid)

    GC.@preserve _core_ptr begin
        @ccall libbraid.braid_Init(
            _app.comm::MPI.MPI_Comm, _app.comm_t::MPI.MPI_Comm,
            tstart::Cdouble, tstop::Cdouble, ntime::Cint,
            _app::Ref{BraidApp}, _c_step::Ptr{Cvoid}, _c_init::Ptr{Cvoid}, _c_clone::Ptr{Cvoid},
            _c_free::Ptr{Cvoid}, _c_sum::Ptr{Cvoid}, _c_norm::Ptr{Cvoid}, _c_access::Ptr{Cvoid},
            _c_bufsize::Ptr{Cvoid}, _c_bufpack::Ptr{Cvoid}, _c_bufunpack::Ptr{Cvoid},
            _core_ptr::Ptr{Ptr{Cvoid}},
        )::Cint
    end

    _core = unsafe_load(_core_ptr)
    Base.Libc.free(_core_ptr)

    if sync !== nothing
        GC.@preserve _core begin
            @ccall libbraid.braid_SetSync(
                _core::Ptr{Cvoid}, _c_sync::Ptr{Cvoid}
            )::Cint
        end
    end

    return BraidCore(_core, _app, tstart, tstop, ntime)
end

"""
Warmup the BraidCore object. XBraid.Drive calls this function automatically, but it can be called manually to precompile all user functions.

This function is not necessary for XBraid to run, but it can significantly reduce the time to run the first time step. Note that this function calls all user functions, so if any user functions have side-effects, this may behave unexpectedly

julia> XBraid.Warmup(core)

Julia> XBraid.Drive(core; warmup=false)

See also: XBraid.Drive
"""
function Warmup(core::BraidCore)
    _app = core._braid_app
    # precompile all user functions by calling them from braid_Warmup
    # if any user functions have side-effects, this may behave unexpectedly
    fdt = (core.tstop - core.tstart) / (core.ntime + 1)
    cdt = 2fdt

    if (_app.basis_init !== nothing) && (_app.inner_prod !== nothing)
        @ccall libbraid.braid_Warmup(
            _app::Ref{BraidApp}, _app.comm::MPI.MPI_Comm, core.tstart::Cdouble, fdt::Cdouble, cdt::Cdouble,
            _c_init::Ptr{Cvoid}, _c_access::Ptr{Cvoid}, _c_free::Ptr{Cvoid},
            _c_clone::Ptr{Cvoid}, _c_sum::Ptr{Cvoid}, _c_norm::Ptr{Cvoid},
            _c_bufsize::Ptr{Cvoid}, _c_bufpack::Ptr{Cvoid}, _c_bufunpack::Ptr{Cvoid},
            C_NULL::Ptr{Cvoid}, C_NULL::Ptr{Cvoid}, _c_step::Ptr{Cvoid},
            _c_init_basis::Ptr{Cvoid}, _c_inner_prod::Ptr{Cvoid}
        )::Cint
    else
        @ccall libbraid.braid_Warmup(
            _app::Ref{BraidApp}, _app.comm::MPI.MPI_Comm, core.tstart::Cdouble, fdt::Cdouble, cdt::Cdouble,
            _c_init::Ptr{Cvoid}, _c_access::Ptr{Cvoid}, _c_free::Ptr{Cvoid},
            _c_clone::Ptr{Cvoid}, _c_sum::Ptr{Cvoid}, _c_norm::Ptr{Cvoid},
            _c_bufsize::Ptr{Cvoid}, _c_bufpack::Ptr{Cvoid}, _c_bufunpack::Ptr{Cvoid},
            C_NULL::Ptr{Cvoid}, C_NULL::Ptr{Cvoid}, _c_step::Ptr{Cvoid},
            C_NULL::Ptr{Cvoid}, C_NULL::Ptr{Cvoid}
        )::Cint
    end
    nothing
end

"""
Wraps the XBraid braid_Drive function. This function is called by XBraid.Drive, and should not be called directly.
"""
function _Drive(core::BraidCore)
    GC.@preserve core begin
        @ccall libbraid.braid_Drive(core._braid_core::Ptr{Cvoid})::Cint
    end
end

"""
Run the XBraid solver. This function calls XBraid.Warmup by default, but this can be disabled by setting warmup=false, in which case a less expensive (but less effective) option is used to precompile the user functions.
"""
function Drive(core::BraidCore; warmup=true)
    if warmup
        Warmup(core)
    else
        _app = core._braid_app
        # cheaper precompile option, but less effective
        begin
            postInitPrecompile(_app)
            if (_app.basis_init !== nothing) && (_app.inner_prod !== nothing)
                deltaPrecompile(_app)
            end
        end
    end
    _Drive(core)
end

function printStats(core::BraidCore)
    GC.@preserve core begin
        @ccall libbraid.braid_PrintStats(core._braid_core::Ptr{Cvoid})::Cint
    end
end

function setTimerFile(core::BraidCore, filestem::String)
    len = @ccall strlen(filestem::Cstring)::Cint
    GC.@preserve core begin
        @ccall libbraid.braid_SetTimerFile(core._braid_core::Ptr{Cvoid}, len::Cint, filestem::Cstring)::Cint
    end
end

function printTimers(core::BraidCore)
    GC.@preserve core begin
        @ccall libbraid.braid_PrintTimers(core._braid_core::Ptr{Cvoid})::Cint
    end
end

function resetTimer(core::BraidCore)
    GC.@preserve core begin
        @ccall libbraid.braid_ResetTimer(core._braid_core::Ptr{Cvoid})::Cint
    end
end

function setTimings(core::BraidCore, timing_level::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetTimings(core._braid_core::Ptr{Cvoid}, timing_level::Cint)::Cint
    end
end

function writeConvHistory(core::BraidCore, filename::String)
    GC.@preserve core begin
        @ccall libbraid.braid_WriteConvHistory(core._braid_core::Ptr{Cvoid}, filename::Cstring)::Cint
    end
end

function setMaxLevels(core::BraidCore, max_levels::Integer)
    @assert max_levels > 0 "Max. levels must be an integer greater than 0"
    GC.@preserve core begin
        @ccall libbraid.braid_SetMaxLevels(core._braid_core::Ptr{Cvoid}, max_levels::Cint)::Cint
    end
end

function setIncrMaxLevels(core::BraidCore)
    GC.@preserve core begin
        @ccall libbraid.braid_SetIncrMaxLevels(core._braid_core::Ptr{Cvoid})::Cint
    end
end

function setSkip(core::BraidCore, skip::Bool)
    GC.@preserve core begin
        @ccall libbraid.braid_SetSkip(core._braid_core::Ptr{Cvoid}, skip::Cint)::Cint
    end
end

function setRefine(core::BraidCore, refine::Bool)
    GC.@preserve core begin
        @ccall libbraid.braid_SetRefine(core._braid_core::Ptr{Cvoid}, refine::Cint)::Cint
    end
end

function setMaxRefinements(core::BraidCore, max_refinements::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetMaxRefinements(core._braid_core::Ptr{Cvoid}, max_refinements::Cint)::Cint
    end
end

function setTPointsCutoff(core::BraidCore, tpoints_cutoff::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetTPointsCutoff(core._braid_core::Ptr{Cvoid}, tpoints_cutoff::Cint)::Cint
    end
end


function setMinCoarse(core::BraidCore, min_coarse::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetMinCoarse(core._braid_core::Ptr{Cvoid}, min_coarse::Cint)::Cint
    end
end

function setRelaxOnlyCG(core::BraidCore, relax_only_cg::Bool)
    GC.@preserve core begin
        @ccall libbraid.braid_SetRelaxOnlyCG(core._braid_core::Ptr{Cvoid}, relax_only_cg::Cint)::Cint
    end
end

function setAbsTol(core::BraidCore, atol::Real)
    GC.@preserve core begin
        @ccall libbraid.braid_SetAbsTol(core._braid_core::Ptr{Cvoid}, atol::Cdouble)::Cint
    end
end

function setRelTol(core::BraidCore, rtol::Real)
    GC.@preserve core begin
        @ccall libbraid.braid_SetRelTol(core._braid_core::Ptr{Cvoid}, rtol::Cdouble)::Cint
    end
end

function setNRelax(core::BraidCore, level::Integer, nrelax::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetNRelax(core._braid_core::Ptr{Cvoid}, level::Cint, nrelax::Cint)::Cint
    end
end

function setCRelaxWt(core::BraidCore, level::Integer, Cwt::Real)
    GC.@preserve core begin
        @ccall libbraid.braid_SetCRelaxWt(core._braid_core::Ptr{Cvoid}, level::Cint, Cwt::Cdouble)::Cint
    end
end

function setCFactor(core::BraidCore, level::Integer, cfactor::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetCFactor(core._braid_core::Ptr{Cvoid}, level::Cint, cfactor::Cint)::Cint
    end
end

function setMaxIter(core::BraidCore, max_iter::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetMaxIter(core._braid_core::Ptr{Cvoid}, max_iter::Cint)::Cint
    end
end

function setFMG(core::BraidCore)
    GC.@preserve core begin
        @ccall libbraid.braid_SetFMG(core._braid_core::Ptr{Cvoid})::Cint
    end
end

function setNFMG(core::BraidCore, k::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetNFMG(core._braid_core::Ptr{Cvoid}, k::Cint)::Cint
    end
end

function setNFMGVcyc(core::BraidCore, nfmg_Vcyc::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetNFMGVcyc(core._braid_core::Ptr{Cvoid}, nfmg_Vcyc::Cint)::Cint
    end
end

function setStorage(core::BraidCore, storage::Bool)
    GC.@preserve core begin
        @ccall libbraid.braid_SetStorage(core._braid_core::Ptr{Cvoid}, storage::Cint)::Cint
    end
end

function setTemporalNorm(core::BraidCore, tnorm::Bool)
    GC.@preserve core begin
        @ccall libbraid.braid_SetTemporalNorm(core._braid_core::Ptr{Cvoid}, tnorm::Cint)::Cint
    end
end

function setPeriodic(core::BraidCore, periodic::Bool)
    GC.@preserve core begin
        @ccall libbraid.braid_SetPeriodic(core._braid_core::Ptr{Cvoid}, periodic::Cint)::Cint
    end
end

function setPrintLevel(core::BraidCore, print_level::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetPrintLevel(core._braid_core::Ptr{Cvoid}, print_level::Cint)::Cint
    end
end

function setFileIOLevel(core::BraidCore, io_level::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetFileIOLevel(core._braid_core::Ptr{Cvoid}, io_level::Cint)::Cint
    end
end

function setPrintFile(core::BraidCore, filename::String)
    GC.@preserve core begin
        @ccall libbraid.braid_SetPrintFile(core._braid_core::Ptr{Cvoid}, filename::Cstring)::Cint
    end
end

function setDefaultPrintFile(core::BraidCore)
    GC.@preserve core begin
        @ccall libbraid.braid_SetDefaultPrintFile(core._braid_core::Ptr{Cvoid})::Cint
    end
end

"""
Set access level for XBraid. This controls how often the user's access function is called.

	- 0: Never call access function
	- 1: Call access function only when XBraid is finished
	- 2: Call access function at every XBraid iteration, on every level

Default is 1.
"""
function setAccessLevel(core::BraidCore, access_level::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetAccessLevel(core._braid_core::Ptr{Cvoid}, access_level::Cint)::Cint
    end
end

function setFinalFCRelax(core::BraidCore)
    GC.@preserve core begin
        @ccall libbraid.braid_SetFinalFCRelax(core._braid_core::Ptr{Cvoid})::Cint
    end
end

function getNumIter(core::BraidCore)
    niter = Ref{Cint}(0)
    GC.@preserve core begin
        @ccall libbraid.braid_GetNumIter(core._braid_core::Ptr{Cvoid}, niter::Ref{Cint})::Cint
    end
    return niter[]
end

function getRNorms(core::BraidCore)
    nrequest = Ref{Cint}()
    nrequest[] = getNumIter(core)
    rnorms = zeros(nrequest[])
    GC.@preserve core begin
        @ccall libbraid.braid_GetRNorms(core._braid_core::Ptr{Cvoid}, nrequest::Ref{Cint}, rnorms::Ref{Cdouble})::Cint
    end
    nrequest[] == 0 && return Float64[]
    return rnorms[1:nrequest[]]
end

function getNLevels(core::BraidCore)
    nlevels = Ref(0)
    GC.@preserve core begin
        @ccall libbraid.braid_GetNLevels(core._braid_core::Ptr{Cvoid}, nlevels::Ref{Cint})::Cint
    end
    return nlevels[]
end

function setSeqSoln(core::BraidCore, seq_soln::Bool)
    GC.@preserve core begin
        @ccall libbraid.braid_GetNLevels(core._braid_core::Ptr{Cvoid}, seq_soln::Cint)::Cint
    end
end

function setRichardsonEstimation(core::BraidCore, est_error::Bool, richardson::Bool, local_order::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetRichardsonEstimation(core._braid_core::Ptr{Cvoid}, est_error::Cint, richardson::Cint, local_order::Cint)::Cint
    end
end

function setSync(core::BraidCore, sync::Function)
    core._braid_app.sync = sync

    GC.@preserve core begin
        @ccall libbraid.braid_SetSync(core._braid_core::Ptr{Cvoid}, _c_sync::Ptr{Cvoid})::Cint
    end
end

function setDeltaCorrection(core::BraidCore, rank::Integer, basis_init::Function; inner_prod::Function=default_inner_prod)
    core._braid_app.basis_init = basis_init
    core._braid_app.inner_prod = inner_prod
    # deltaPrecompile(core._braid_app)

    GC.@preserve core begin
        @ccall libbraid.braid_SetDeltaCorrection(core._braid_core::Ptr{Cvoid}, rank::Cint, _c_init_basis::Ptr{Cvoid}, _c_inner_prod::Ptr{Cvoid})::Cint
    end
end

function setDeferDelta(core::BraidCore, level::Integer, iter::Integer)
    GC.@preserve core begin
        @ccall libbraid.braid_SetDeferDelta(core._braid_core::Ptr{Cvoid}, level::Cint, iter::Cint)::Cint
    end
end

function setLyapunovEstimation(core::BraidCore; relax::Bool=false, cglv::Bool=true, exponents::Bool=false)
    GC.@preserve core begin
        @ccall libbraid.braid_SetLyapunovEstimation(core._braid_core::Ptr{Cvoid}, relax::Cint, cglv::Cint, exponents::Cint)::Cint
    end
end

# Still missing spatial coarsening, residual, and adjoint

# braid_test
function testInitAccess(app::BraidApp, t::Real, outputFile::Libc.FILE)
    pass = @ccall libbraid.braid_TestInitAccess(
        app::Ref{BraidApp},
        app.comm::MPI.MPI_Comm,
        outputFile::Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_access::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
    )::Cint
    return pass > 0

    # println("Serialized size of user vector: $(app.bufsize)")
    # println("Check output for objects not properly freed:")
    # println(app.ref_ids)
end
testInitAccess(app::BraidApp, t::Real, outputFile::IO) = testInitAccess(app, t, Libc.FILE(outputFile))
testInitAccess(app::BraidApp, t::Real) = testInitAccess(app, t, c_stdout)

function testClone(app::BraidApp, t::Real, outputFile::Libc.FILE)
    pass = @ccall libbraid.braid_TestClone(
        app::Ref{BraidApp},
        app.comm::MPI.MPI_Comm,
        outputFile::Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_access::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
        _c_clone::Ptr{Cvoid},
    )::Cint
    return pass > 0
end
testClone(app::BraidApp, t::Real, outputFile::IO) = testClone(app, t, Libc.FILE(outputFile))
testClone(app::BraidApp, t::Real) = testClone(app, t, c_stdout)

function testSum(app::BraidApp, t::Real, outputFile::Libc.FILE)
    pass = @ccall libbraid.braid_TestSum(
        app::Ref{BraidApp},
        app.comm::MPI.MPI_Comm,
        outputFile::Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_access::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
        _c_clone::Ptr{Cvoid},
        _c_sum::Ptr{Cvoid},
    )::Cint
    return pass > 0
end
testSum(app::BraidApp, t::Real, outputFile::IO) = testSum(app, t, Libc.FILE(outputFile))
testSum(app::BraidApp, t::Real) = testSum(app, t, c_stdout)

function testSpatialNorm(app::BraidApp, t::Real, outputFile::Libc.FILE)
    pass = @ccall libbraid.braid_TestSpatialNorm(
        app::Ref{BraidApp},
        app.comm::MPI.MPI_Comm,
        outputFile::Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
        _c_clone::Ptr{Cvoid},
        _c_sum::Ptr{Cvoid},
        _c_norm::Ptr{Cvoid},
    )::Cint
    return pass > 0
end
testSpatialNorm(app::BraidApp, t::Real, outputFile::IO) = testSpatialNorm(app, t, Libc.FILE(outputFile))
testSpatialNorm(app::BraidApp, t::Real) = testSpatialNorm(app, t, c_stdout)

function testBuf(app::BraidApp, t::Real, outputFile::Libc.FILE)
    pass = @ccall libbraid.braid_TestBuf(
        app::Ref{BraidApp},
        app.comm::MPI.MPI_Comm,
        outputFile::Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
        _c_sum::Ptr{Cvoid},
        _c_norm::Ptr{Cvoid},
        _c_bufsize::Ptr{Cvoid},
        _c_bufpack::Ptr{Cvoid},
        _c_bufunpack::Ptr{Cvoid},
    )::Cint
    return pass > 0
    # print('\n')
    # println("Check output for objects not freed:")
    # println(app.ref_ids)
end
testBuf(app::BraidApp, t::Real, outputFile::IO) = testBuf(app, t, Libc.FILE(outputFile))
testBuf(app::BraidApp, t::Real) = testBuf(app, t, c_stdout)

function testDelta(app::BraidApp, t::Real, dt::Real, rank::Integer, outputFile::Libc.FILE)
    app.basis_init::Function
    app.inner_prod::Function
    pass = @ccall libbraid.braid_TestDelta(
        app::Ref{BraidApp},
        app.comm::MPI.MPI_Comm,
        outputFile::Libc.FILE,
        t::Cdouble,
        dt::Cdouble,
        rank::Cint,
        _c_init::Ptr{Cvoid},
        _c_init_basis::Ptr{Cvoid},
        _c_access::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
        _c_clone::Ptr{Cvoid},
        _c_sum::Ptr{Cvoid},
        _c_bufsize::Ptr{Cvoid},
        _c_bufpack::Ptr{Cvoid},
        _c_bufunpack::Ptr{Cvoid},
        _c_inner_prod::Ptr{Cvoid},
        _c_step::Ptr{Cvoid},
    )::Cint
    return pass > 0
    # println("Check output for objects not freed:")
    # println(app.ref_ids)
end
testDelta(app::BraidApp, t::Real, dt::Real, rank::Integer, outputFile::IO) = testDelta(app, t, dt, rank, Libc.FILE(outputFile))
testDelta(app::BraidApp, t::Real, dt::Real, rank::Integer) = testDelta(app, t, dt, rank, c_stdout)

end # module XBraid
