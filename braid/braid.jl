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

using Serialization: serialize, deserialize # used to pack arbitrary julia objects into buffers
using MPI
C_stdout = Libc.FILE(Libc.RawFD(1), "w")  # corresponds to C standard output
XBraid = "./libbraid.so"

#=
 |  For this to work, XBraid must be compiled as a shared library, using 
 |  $ cd braid/
 |  $ make libbraid.so
 |  and the Fortran flags braid_Fortran_SpatialCoarsen, braid_Fortran_Residual,
 |  braid_Fortran_TimeGrid, and braid_Fortran_Sync must all be set to zero in
 |  braid.h before compiling. (TODO: figure out why this is...)
 =#

#= 
 | This interface uses the braid_app as an internal structure that is not exposed to the user.
 | This is so that the user interface only uses memory safe function calls
 | User defined data structures that are not time-dependent should be stored in globally scoped variables
 =#
mutable struct _jl_braid_app
    comm::MPI.Comm    # global mpi communicator
    comm_t::MPI.Comm  # temporal mpi communicator

    # required user functions
    step::Function
    init::Function
    sum::Function
    spatialnorm::Function
    access::Function

    # optional user functions
    basis_init::Union{Function, Nothing}
    inner_prod::Union{Function, Nothing}

    ref_ids::IdDict   # dictionary to store globally scoped references to allocated braid_vector objects
    bufsize::Integer  # expected serialized size of user's braid_vector object
    bufsize_lyap::Integer
end


mutable struct _jl_braid_core
end

#= 
 | stores internal structures _braid_core and _braid_app,
 =#
struct braid_core
    # internal values
    _braid_core::Ptr{_jl_braid_core}
    _braid_app::_jl_braid_app

end

# This can contain anything
mutable struct _jl_braid_vector
    user_vector
end

# some constructors
function _jl_braid_app(comm::MPI.Comm, step::Function, init::Function, sum::Function, norm::Function, access::Function) 
    _jl_braid_app(comm, comm, step, init, sum, norm, access, nothing, nothing, IdDict(), 0, 0)
end

function _jl_braid_app(comm::MPI.Comm, comm_t::MPI.Comm, step::Function, init::Function, sum::Function, norm::Function, access::Function) 
    _jl_braid_app(comm, comm_t, step, init, sum, norm, access, nothing, nothing, IdDict(), 0, 0)
end

# Used to add/remove a reference to the vector u to the IdDict stored in the app.
# this keeps all newly allocated braid_vectors in the global scope, preventing them
# from being garbage collected!
function _register_vector(app::_jl_braid_app, u::_jl_braid_vector)
    app.ref_ids[objectid(u)] = u
end

function _deregister_vector(app::_jl_braid_app, u::_jl_braid_vector)
    pop!(app.ref_ids, objectid(u))
end

# these are internal functions which directly interface with XBraid
# TODO: get the status structures working
function _jl_step!(_app::Ptr{Cvoid},
                   _ustop::Ptr{Cvoid},
                   _fstop::Ptr{Cvoid},
                   _u::Ptr{Cvoid},
                   status::Ptr{Cvoid})::Cint
    # println("step")
    app = unsafe_pointer_to_objref(_app)
    u = unsafe_pointer_to_objref(_u)
    ustop = unsafe_pointer_to_objref(_ustop)
    tstart, tstop = Ref(0.0), Ref(0.0) # guaranteed to not be garbage collected until dereferenced
    @ccall XBraid.braid_StepStatusGetTstartTstop(status::Ptr{Cvoid}, tstart::Ref{Cdouble}, tstop::Ref{Cdouble})::Cint
    if _fstop !== C_NULL
        fstop = unsafe_pointer_to_objref(_fstop)
        app.step(u.user_vector, ustop.user_vector, fstop.user_vector, tstart[], tstop[])
    else
        app.step(u.user_vector, ustop.user_vector, nothing, tstart[], tstop[])
    end

    return 0
end
_c_step = @cfunction(_jl_step!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))


function _jl_init!(_app::Ptr{Cvoid}, t::Cdouble, u_ptr::Ptr{Ptr{Cvoid}})::Cint
    # println("init")
    app = unsafe_pointer_to_objref(_app)
    # initialize u and register a reference with IdDict
    u = _jl_braid_vector(app.init(t))
    _register_vector(app, u)

    unsafe_store!(u_ptr, pointer_from_objref(u))

    # store max size of all initialized vectors
    if app.bufsize == 0
        buffer = IOBuffer()
        serialize(buffer, u)
        app.bufsize = buffer.ptr
    end

    # TODO: figure out how to actually do this properly
    # u_size = Base.summarysize(Ref(u)) + 9
    # if u_size > app.bufsize
    #     app.bufsize = u_size
    # end

    return 0
end
_c_init = @cfunction(_jl_init!, Cint, (Ptr{Cvoid}, Cdouble, Ptr{Ptr{Cvoid}}))

function _jl_init_basis!(_app::Ptr{Cvoid}, t::Cdouble, index::Cint, u_ptr::Ptr{Ptr{Cvoid}})::Cint
    app = unsafe_pointer_to_objref(_app)
    u = _jl_braid_vector(app.init_basis(t, index))
    _register_vector(app, u)
    unsafe_store!(u_ptr, pointer_from_objref(u))

    # store max size of all initialized vectors
    if app.bufsize == 0
        buffer = IOBuffer()
        serialize(buffer, u)
        app.bufsize_lyap = buffer.ptr
    end

    return 0
end
_c_init_basis = @cfunction(_jl_init_basis!, Cint, (Ptr{Cvoid}, Cdouble, Cint, Ptr{Ptr{Cvoid}}))

function _jl_clone!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, v_ptr::Ptr{Ptr{Cvoid}})::Cint
    # println("clone")
    app = unsafe_pointer_to_objref(_app)
    u = unsafe_pointer_to_objref(_u)
    # initialize v, and copy u into v
    v = deepcopy(u)

    # then register v with IdDict and store in v_ptr
    _register_vector(app, v)
    unsafe_store!(v_ptr, pointer_from_objref(v))

    return 0
end
_c_clone = @cfunction(_jl_clone!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}))

function _jl_free!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid})::Cint
    # println("free")
    app = unsafe_pointer_to_objref(_app)
    u = unsafe_pointer_to_objref(_u)
    # removing the global reference to u will cause it to be garbage collected
    _deregister_vector(app, u)
    u = C_NULL
    return 0
end
_c_free = @cfunction(_jl_free!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))

function _jl_sum!(_app::Ptr{Cvoid},
    alpha::Cdouble, _x::Ptr{Cvoid},
    beta::Cdouble, _y::Ptr{Cvoid})::Cint
    # println("sum")
    app = unsafe_pointer_to_objref(_app)
    x = unsafe_pointer_to_objref(_x)
    y = unsafe_pointer_to_objref(_y)
    app.sum(alpha, x.user_vector, beta, y.user_vector)
    # y.user_vector = alpha * x.user_vector + beta * y.user_vector
    return 0
end
_c_sum = @cfunction(_jl_sum!, Cint, (Ptr{Cvoid}, Cdouble, Ptr{Cvoid}, Cdouble, Ptr{Cvoid}))

function _jl_norm!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, norm_ptr::Ptr{Cdouble})::Cint
    # println("norm")
    app = unsafe_pointer_to_objref(_app)
    u = unsafe_pointer_to_objref(_u)
    norm = app.spatialnorm(u.user_vector)
    unsafe_store!(norm_ptr, norm)

    return 0
end
_c_norm = @cfunction(_jl_norm!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))

function _jl_inner_prod!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, _v::Ptr{Cvoid}, norm_ptr::Ptr{Cdouble})::Cint
    app = unsafe_pointer_to_objref(_app)
    u = unsafe_pointer_to_objref(_u)
    v = unsafe_pointer_to_objref(_v)
    # TODO: call user's norm function
    prod = app.inner_prod(u.user_vector, v.user_vector)
    unsafe_store!(norm_ptr, prod)

    return 0
end
_c_inner_prod = @cfunction(_jl_inner_prod!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))

function _jl_access!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, status::Ptr{Cvoid})::Cint
    # println("access")
    app = unsafe_pointer_to_objref(_app)
    u = unsafe_pointer_to_objref(_u)
    if !isnothing(app.access)
        app.access(u.user_vector)
    end
    return 0
end
_c_access = @cfunction(_jl_access!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))

function _jl_bufsize!(_app::Ptr{Cvoid}, size_ptr::Ptr{Cint}, status::Ptr{Cvoid})::Cint
    # println("bufsize")
    app = unsafe_pointer_to_objref(_app)
    unsafe_store!(size_ptr, app.bufsize)
    return 0
end
_c_bufsize = @cfunction(_jl_bufsize!, Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cvoid}))

function _jl_bufpack!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, _buffer::Ptr{Cvoid}, status::Ptr{Cvoid})::Cint
    # println("bufpack")
    app = unsafe_pointer_to_objref(_app)
    u = unsafe_pointer_to_objref(_u)
    buff_arr = unsafe_wrap(Vector{UInt8}, Base.unsafe_convert(Ptr{UInt8}, _buffer), app.bufsize)
    buffer = IOBuffer(buff_arr, write=true, maxsize=app.bufsize)
    serialize(buffer, u)
    # println(buffer.data)
    # show(buffer)
    # print('\n')

    # buffer = IOBuffer()
    # serialize(buffer, u)
    # println(buffer.data)
    # show(buffer)
    # print('\n')

    # tell XBraid the written size
    @ccall XBraid.braid_BufferStatusSetSize(status::Ptr{Cvoid}, buffer.size::Cdouble)::Cint

    return 0
end
_c_bufpack = @cfunction(_jl_bufpack!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))

function _jl_bufunpack!(_app::Ptr{Cvoid}, _buffer::Ptr{Cvoid}, u_ptr::Ptr{Ptr{Cvoid}}, status::Ptr{Cvoid})::Cint
    # println("bufunpack")
    app = unsafe_pointer_to_objref(_app)
    buff_arr = unsafe_wrap(Vector{UInt8}, Base.unsafe_convert(Ptr{UInt8}, _buffer), app.bufsize)
    buffer = IOBuffer(buff_arr, read=true, write=true, maxsize=app.bufsize)
    # println(buffer.data)
    # show(buffer)
    # print('\n')

    # unpack the buffer into a new julia object, then register with IdDict
    u = deserialize(buffer)
    _register_vector(app, u)
    # store u in provided pointer
    unsafe_store!(u_ptr, pointer_from_objref(u))
    return 0
end
_c_bufunpack = @cfunction(_jl_bufunpack!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))

function Init(comm_world::MPI.Comm, comm_t::MPI.Comm,
              tstart::Real, tstop::Real, ntime::Int, step::Function,
              init::Function, sum::Function, spatialnorm::Function, access::Function)::braid_core
    _app = _jl_braid_app(comm_world, comm_t, step, init, sum, spatialnorm, access)

    _core_ptr = reinterpret(Ptr{Ptr{_jl_braid_core}}, pointer_from_objref(Ref(C_NULL)))
    @ccall XBraid.braid_Init(_app.comm::MPI.MPI_Comm, _app.comm_t::MPI.MPI_Comm, tstart::Cdouble, tstop::Cdouble, ntime::Cint,
                             _app::Ref{_jl_braid_app}, _c_step::Ptr{Cvoid}, _c_init::Ptr{Cvoid}, _c_clone::Ptr{Cvoid},
                             _c_free::Ptr{Cvoid}, _c_sum::Ptr{Cvoid}, _c_norm::Ptr{Cvoid}, _c_access::Ptr{Cvoid},
                             _c_bufsize::Ptr{Cvoid}, _c_bufpack::Ptr{Cvoid}, _c_bufunpack::Ptr{Cvoid},
                             _core_ptr::Ptr{Ptr{_jl_braid_core}})::Cint

    _core = unsafe_load(_core_ptr)
    return braid_core(_core, _app)
end

function Drive(core::braid_core)
    @ccall XBraid.braid_Drive(core._braid_core::Ptr{_jl_braid_core})::Cint
end

# this lets XBraid deallocate anything that was allocated in C
function Destroy!(core::braid_core)
    @ccall XBraid.braid_Destroy(core._braid_core::Ptr{_jl_braid_core})::Cint
    core = nothing # julia will destroy the rest
end

function PrintStats(core::braid_core)
    @ccall XBraid.braid_PrintStats(core._braid_core::Ptr{_jl_braid_core})::Cint
end

function SetTimerFile(core::braid_core, filestem::String)
    len = @ccall strlen(filestem::Cstring)::Cint
    @ccall XBraid.braid_SetTimerFile(core._braid_core::Ptr{_jl_braid_core}, len::Cint, filestem::Cstring)::Cint
end

PrintTimers(core::braid_core) = @ccall XBraid.braid_SetTimerFile(core._braid_core::Ptr{_jl_braid_core})::Cint
ResetTimer(core::braid_core) = @ccall XBraid.braid_ResetTimer(core._braid_core::Ptr{_jl_braid_core})::Cint
braid_SetTimings(core::braid_core, timing_level::Integer) = @ccall XBraid.braid_SetTimings(core._braid_core::Ptr{_jl_braid_core}, timing_level::Cint)::Cint

function WriteConvHistory(core::braid_core, filename::String)
    @ccall XBraid.braid_WriteConvHistory(core._braid_core::Ptr{_jl_braid_core}, filename::Cstring)::Cint
end

function SetMaxLevels(core::braid_core, max_levels::Integer)
    @ccall XBraid.braid_SetMaxLevels(core._braid_core::Ptr{_jl_braid_core}, max_levels::Cint)::Cint
end

SetIncrMaxLevels(core::braid_core) = @ccall XBraid.braid_SetIncrMaxLevels(core._braid_core::Ptr{_jl_braid_core})::Cint

SetSkip(core::braid_core, skip::Bool) = @ccall XBraid.braid_SetSkip(core._braid_core::Ptr{_jl_braid_core}, skip::Cint)::Cint
SetRefine(core::braid_core, refine::Bool) = @ccall XBraid.braid_SetSkip(core._braid_core::Ptr{_jl_braid_core}, refine::Cint)::Cint
SetMaxRefinements(core::braid_core, max_refinements::Integer) = @ccall XBraid.braid_SetMaxRefinements(core._braid_core::Ptr{_jl_braid_core}, max_refinements::Cint)::Cint
SetTPointsCutoff(core::braid_core, tpoints_cutoff::Integer) = @ccall XBraid.braid_SetTPointsCutoff(core._braid_core::Ptr{_jl_braid_core}, tpoints_cutoff::Cint)::Cint
SetMinCoarse(core::braid_core, min_coarse::Integer) = @ccall XBraid.braid_SetMinCoarse(core._braid_core::Ptr{_jl_braid_core}, min_coarse::Cint)::Cint
SetRelaxOnlyCG(core::braid_core, relax_only_cg::Bool) = @ccall XBraid.braid_SetRelaxOnlyCG(core._braid_core::Ptr{_jl_braid_core}, relax_only_cg::Cint)::Cint
SetAbsTol(core::braid_core, atol::Real) = @ccall XBraid.braid_SetAbsTol(core._braid_core::Ptr{_jl_braid_core}, atol::Cdouble)::Cint
SetRelTol(core::braid_core, rtol::Real) = @ccall XBraid.braid_SetRelTol(core._braid_core::Ptr{_jl_braid_core}, rtol::Cdouble)::Cint
SetNRelax(core::braid_core, level::Integer, nrelax::Integer) = @ccall XBraid.braid_SetNRelax(core._braid_core::Ptr{_jl_braid_core}, level::Cint, nrelax::Cint)::Cint
SetCRelaxWt(core::braid_core, level::Integer, Cwt::Real) = @ccall XBraid.braid_SetCRelaxWt(core._braid_core::Ptr{_jl_braid_core}, level::Cint, Cwt::Cdouble)::Cint
SetCFactor(core::braid_core, level::Integer, cfactor::Integer) = @ccall XBraid.braid_SetNRelax(core._braid_core::Ptr{_jl_braid_core}, level::Cint, cfactor::Cint)::Cint
SetMaxIter(core::braid_core, max_iter::Integer) = @ccall XBraid.braid_SetMaxIter(core._braid_core::Ptr{_jl_braid_core}, max_iter::Cint)::Cint
SetFMG(core::braid_core) = @ccall XBraid.braid_SetFMG(core._braid_core::Ptr{_jl_braid_core})::Cint
SetNFMG(core::braid_core, k::Integer) = @ccall XBraid.braid_SetNFMG(core._braid_core::Ptr{_jl_braid_core}, k::Cint)::Cint
SetNFMGVcyc(core::braid_core, nfmg_Vcyc::Integer) = @ccall XBraid.braid_SetNFMGVcyc(core._braid_core::Ptr{_jl_braid_core}, nfmg_Vcyc::Cint)::Cint
SetStorage(core::braid_core, storage::Bool) = @ccall XBraid.braid_SetStorage(core._braid_core::Ptr{_jl_braid_core}, storage::Cint)::Cint
SetTemporalNorm(core::braid_core, tnorm::Bool) = @ccall XBraid.braid_SetTemporalNorm(core._braid_core::Ptr{_jl_braid_core}, tnorm::Cint)::Cint
SetPeriodic(core::braid_core, periodic::Bool) = @ccall XBraid.braid_SetPeriodic(core._braid_core::Ptr{_jl_braid_core}, periodic::Cint)::Cint
SetPrintLevel(core::braid_core, print_level::Integer) = @ccall XBraid.braid_SetPrintLevel(core._braid_core::Ptr{_jl_braid_core}, print_level::Cint)::Cint
SetFileIOLevel(core::braid_core, io_level::Integer) = @ccall XBraid.braid_SetFileIOLevel(core._braid_core::Ptr{_jl_braid_core}, io_level::Cint)::Cint
SetPrintFile(core::braid_core, filename::String) = @ccall XBraid.braid_SetPrintFile(core._braid_core::Ptr{_jl_braid_core}, filename::Cstring)::Cint
SetDefaultPrintFile(core::braid_core) = @ccall XBraid.braid_SetDefaultPrintFile(core._braid_core::Ptr{_jl_braid_core})::Cint
SetAccessLevel(core::braid_core, access_level::Integer) = @ccall XBraid.braid_SetAccessLevel(core._braid_core::Ptr{_jl_braid_core}, access_level::Cint)::Cint
SetFinalFCRelax(core::braid_core) = @ccall XBraid.braid_SetFinalFCRelax(core._braid_core::Ptr{_jl_braid_core})::Cint
function GetNumIter(core::braid_core)
    niter = Ref(0)
    @ccall XBraid.braid_GetNumIter(core._braid_core::Ptr{_jl_braid_core}, niter::Ref{Cint})::Cint
    return niter[] # dereference
end

function GetRNorms(core::braid_core)
    nrequest = Ref(GetNumIter(core))
    rnorms = zeros(nrequest)
    @ccall XBraid.braid_GetRNorms(core._braid_core::Ptr{_jl_braid_core}, nrequest::Ref{Cint}, rnorms::Ref{Cdouble})::Cint
    return rnorms
end

function GetNLevels(core::braid_core)
    nlevels = Ref(0)
    @ccall XBraid.braid_GetNLevels(core._braid_core::Ptr{_jl_braid_core}, nlevels::Ref{Cint})::Cint
    return nlevels[]
end

SetSeqSoln(core::braid_core, seq_soln::Bool) = @ccall XBraid.braid_GetNLevels(core._braid_core::Ptr{_jl_braid_core}, seq_soln::Cint)::Cint

function SetRichardsonEstimation(core::braid_core, est_error::Bool, richardson::Bool, local_order::Integer)
    @ccall XBraid.braid_SetRichardsonEstimation(core._braid_core::Ptr{_jl_braid_core}, est_error::Cint, richardson::Cint, local_order::Cint)::Cint
end

function SetDeltaCorrection(core::braid_core, rank::Integer, basis_init::Function, inner_prod::Function)
    core._braid_app.basis_init = basis_init
    core._braid_app.inner_prod = inner_prod
    @ccall XBraid.braid_SetDeltaCorrection(core._braid_core::Ptr{_jl_braid_core}, rank::Cint, _c_basis_init::Ptr{Cvoid}, _c_inner_prod::Ptr{Cvoid})::Cint
end

function SetDeferDelta(core::braid_core, level::Integer, iter::Integer)
    @ccall XBraid.braid_SetDeferDelta(core._braid_core::Ptr{_jl_braid_core}, level::Cint, iter::Cint)::Cint
end

function SetLyapunovEstimation(core::braid_core, relax::Bool, cglv::Bool, exponents::Bool)
    @ccall XBraid.braid_SetLyapunovEstimation(core._braid_core::Ptr{_jl_braid_core}, relax::Cint, cglv::Cint, exponents::Bool)::Cint
end

# Still missing spatial coarsening, sync, residual, and adjoint

# braid_test
function TestInitAccess(app::_jl_braid_app, t::Real)
    @ccall XBraid.braid_TestInitAccess(
        app::Ref{_jl_braid_app},
        app.comm::MPI.MPI_Comm,
        C_stdout::Base.Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_access::Ptr{Cvoid},
        _c_free::Ptr{Cvoid}
    )::Cint

    println("Serialized size of user vector: $(app.bufsize)")
    println("Check output for objects not properly freed:")
    println(app.ref_ids)
end

function TestClone(app::_jl_braid_app, t::Real)
    @ccall XBraid.braid_TestClone(
        app::Ref{_jl_braid_app},
        app.comm::MPI.MPI_Comm,
        C_stdout::Base.Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_access::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
        _c_clone::Ptr{Cvoid}
    )::Cint
end

function TestSum(app::_jl_braid_app, t::Real)
    @ccall XBraid.braid_TestSum(
        app::Ref{_jl_braid_app},
        app.comm::MPI.MPI_Comm,
        C_stdout::Base.Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_access::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
        _c_clone::Ptr{Cvoid},
        _c_sum::Ptr{Cvoid}
    )::Cint
end

function TestSpatialNorm(app::_jl_braid_app, t::Real)
    @ccall XBraid.braid_TestSpatialNorm(
        app::Ref{_jl_braid_app},
        app.comm::MPI.MPI_Comm,
        C_stdout::Base.Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
        _c_clone::Ptr{Cvoid},
        _c_sum::Ptr{Cvoid},
        _c_norm::Ptr{Cvoid}
    )::Cint
end

function TestBuf(app::_jl_braid_app, t::Real)
    @ccall XBraid.braid_TestBuf(
        app::Ref{_jl_braid_app},
        app.comm::MPI.MPI_Comm,
        C_stdout::Base.Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
        _c_sum::Ptr{Cvoid},
        _c_norm::Ptr{Cvoid},
        _c_bufsize::Ptr{Cvoid},
        _c_bufpack::Ptr{Cvoid},
        _c_bufunpack::Ptr{Cvoid}
    )::Cint
    print('\n')
    println("Check output for objects not freed:")
    println(app.ref_ids)
end

# main ex-01
MPI.Init()
comm = MPI.COMM_WORLD

function my_step(u, ustop, fstop, tstart, tstop)
    # backward Euler
    u[] = 1. / (1. + tstop-tstart) * u[]
end

function my_init_basis(t, i)
    return Ref(1.)
end

function my_init(t)
    if t == 0.
        u = Ref(1.)
    else
        u = Ref(0.456)
    end
end

function my_sum(a, x, b, y)
    y[] = a*x[] + b*y[]
end

function my_access(u)
    println(u[])
end

my_norm(u) = abs(u[])
my_innerprod(u, v) = u[] * v[]

app = _jl_braid_app(comm, my_step, my_init, my_sum, my_norm, my_access)

TestInitAccess(app, 0.)
TestClone(app, 0.)
TestSum(app, 0.)
TestSpatialNorm(app, 0.)
TestBuf(app, 0.)


ntime = 10
tstart = 0.0
tstop = tstart + ntime / 2.0;
core = Init(comm, comm, tstart, tstop, ntime, my_step, my_init, my_sum, my_norm, my_access)

SetPrintLevel(core, 2)
SetMaxLevels(core, 2)
SetAbsTol(core, 1.e-6)
SetCFactor(core, -1, 2)

Drive(core)

Destroy!(core)

println("I'm amazed you got here!!!")
