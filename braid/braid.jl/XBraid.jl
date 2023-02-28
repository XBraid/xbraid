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
# Not sure I want to export anything...
# export Init, Drive, etc.

# BEGIN MODULE XBraid
using Serialization: serialize, deserialize # used to pack arbitrary julia objects into buffers
using MPI

include("status.jl")
include("wrapper_functions.jl")

#=
 |  For this to work, XBraid must be compiled as a shared library.
 |  $ make braid shared=yes
 |  and the Fortran flags braid_Fortran_SpatialCoarsen, braid_Fortran_Residual,
 |  braid_Fortran_TimeGrid, and braid_Fortran_Sync must all be set to zero in
 |  braid.h before compiling. (TODO: figure out why this is...)
 =#
libbraid = "../libbraid.so"

C_stdout = Libc.FILE(Libc.RawFD(1), "w")  # corresponds to C standard output
function get_null_double_ptr(Type) 
    pp = reinterpret(Ptr{Ptr{Type}}, pointer_from_objref(Ref(0)))
    unsafe_store!(pp, C_NULL)
    return pp
end

"""
This is an internal structure that is not exposed to the user.
This enables the user interface to only use memory safe function calls.
User defined data structures that are not time-dependent can be declared in 
the global scope, or they can be packed into a user defined object and passed to
Init(), in which case the object will be passed as the first argument in every user
defined function call.
_app.user_app is passed as the first argument to every user defined function.
Stores any time-independent data the user may need to compute a time-step.
Large, preallocate arrays should be packed into this object, since operations
involving globally scoped variables are slower than local variables.
struct braid_App end
 """
mutable struct BraidApp
    user_app::Any     # user defined app data structure (can be anything)

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

# default constructor
function BraidApp(app, comm::MPI.Comm, comm_t::MPI.Comm, step::Function, init::Function, sum::Function, norm::Function, access::Function, basis_init::Function, inner_prod::Function)
    BraidApp(app, comm, comm_t, step, init, sum, norm, access, basis_init, inner_prod, IdDict(), 0, 0)
end

function BraidApp(app, comm::MPI.Comm, comm_t::MPI.Comm, step::Function, init::Function, sum::Function, norm::Function, access::Function)
    BraidApp(app, comm, comm_t, step, init, sum, norm, access, nothing, nothing, IdDict(), 0, 0)
end

function BraidApp(app, comm::MPI.Comm, step::Function, init::Function, sum::Function, norm::Function, access::Function)
    BraidApp(app, comm, comm, step, init, sum, norm, access)
end

# This can contain anything (TODO: do I even need this?)
mutable struct BraidVector
    user_vector
end

"""
wraps the pointer to braid status structures,
automatically calling braid_Destroy(core) when
this is garbage collected
"""
mutable struct BraidCore
    # internal values
    _braid_core::Ptr{C_NULL}
    _braid_app::BraidApp
    function BraidCore(_braid_core, _braid_app)
        x = new(_braid_core, _braid_app)
        finalizer(x) do core
            @async println("Destroying BraidCore")
            @ccall libbraid.braid_Destroy(core._braid_core::Ptr{C_NULL})::Cint
        end
    end
end


# Used to add/remove a reference to the vector u to the IdDict stored in the app.
# this keeps all newly allocated braid_vectors in the global scope, preventing them
# from being garbage collected!
function _register_vector(app::BraidApp, u::BraidVector)
    app.ref_ids[objectid(u)] = u
end

function _deregister_vector(app::BraidApp, u::BraidVector)
    pop!(app.ref_ids, objectid(u))
end

function Init(comm_world::MPI.Comm, comm_t::MPI.Comm,
                   tstart::Real, tstop::Real, ntime::Int, 
                   step::Function, init::Function, sum::Function, spatialnorm::Function, access::Function;
                   app=nothing
                   )::BraidCore
    _app = BraidApp(app, comm_world, comm_t, step, init, sum, spatialnorm, access)

    _core_ptr = get_null_double_ptr(C_NULL)
    @ccall libbraid.braid_Init(_app.comm::MPI.MPI_Comm, _app.comm_t::MPI.MPI_Comm, tstart::Cdouble, tstop::Cdouble, ntime::Cint,
                             _app::Ref{BraidApp}, _c_step::Ptr{Cvoid}, _c_init::Ptr{Cvoid}, _c_clone::Ptr{Cvoid},
                             _c_free::Ptr{Cvoid}, _c_sum::Ptr{Cvoid}, _c_norm::Ptr{Cvoid}, _c_access::Ptr{Cvoid},
                             _c_bufsize::Ptr{Cvoid}, _c_bufpack::Ptr{Cvoid}, _c_bufunpack::Ptr{Cvoid},
                             _core_ptr::Ptr{Ptr{C_NULL}})::Cint

    _core = unsafe_load(_core_ptr)
    return BraidCore(_core, _app)
end

function Drive(core::BraidCore)
    @ccall libbraid.braid_Drive(core._braid_core::Ptr{C_NULL})::Cint
end


function PrintStats(core::BraidCore)
    @ccall libbraid.braid_PrintStats(core._braid_core::Ptr{C_NULL})::Cint
end

function SetTimerFile(core::BraidCore, filestem::String)
    len = @ccall strlen(filestem::Cstring)::Cint
    @ccall libbraid.braid_SetTimerFile(core._braid_core::Ptr{C_NULL}, len::Cint, filestem::Cstring)::Cint
end

PrintTimers(core::BraidCore) = @ccall libbraid.braid_SetTimerFile(core._braid_core::Ptr{C_NULL})::Cint
ResetTimer(core::BraidCore) = @ccall libbraid.braid_ResetTimer(core._braid_core::Ptr{C_NULL})::Cint
braid_SetTimings(core::BraidCore, timing_level::Integer) = @ccall libbraid.braid_SetTimings(core._braid_core::Ptr{C_NULL}, timing_level::Cint)::Cint

function WriteConvHistory(core::BraidCore, filename::String)
    @ccall libbraid.braid_WriteConvHistory(core._braid_core::Ptr{C_NULL}, filename::Cstring)::Cint
end

function SetMaxLevels(core::BraidCore, max_levels::Integer)
    @ccall libbraid.braid_SetMaxLevels(core._braid_core::Ptr{C_NULL}, max_levels::Cint)::Cint
end

SetIncrMaxLevels(core::BraidCore) = @ccall libbraid.braid_SetIncrMaxLevels(core._braid_core::Ptr{C_NULL})::Cint

SetSkip(core::BraidCore, skip::Bool) = @ccall libbraid.braid_SetSkip(core._braid_core::Ptr{C_NULL}, skip::Cint)::Cint
SetRefine(core::BraidCore, refine::Bool) = @ccall libbraid.braid_SetSkip(core._braid_core::Ptr{C_NULL}, refine::Cint)::Cint
SetMaxRefinements(core::BraidCore, max_refinements::Integer) = @ccall libbraid.braid_SetMaxRefinements(core._braid_core::Ptr{C_NULL}, max_refinements::Cint)::Cint
SetTPointsCutoff(core::BraidCore, tpoints_cutoff::Integer) = @ccall libbraid.braid_SetTPointsCutoff(core._braid_core::Ptr{C_NULL}, tpoints_cutoff::Cint)::Cint
SetMinCoarse(core::BraidCore, min_coarse::Integer) = @ccall libbraid.braid_SetMinCoarse(core._braid_core::Ptr{C_NULL}, min_coarse::Cint)::Cint
SetRelaxOnlyCG(core::BraidCore, relax_only_cg::Bool) = @ccall libbraid.braid_SetRelaxOnlyCG(core._braid_core::Ptr{C_NULL}, relax_only_cg::Cint)::Cint
SetAbsTol(core::BraidCore, atol::Real) = @ccall libbraid.braid_SetAbsTol(core._braid_core::Ptr{C_NULL}, atol::Cdouble)::Cint
SetRelTol(core::BraidCore, rtol::Real) = @ccall libbraid.braid_SetRelTol(core._braid_core::Ptr{C_NULL}, rtol::Cdouble)::Cint
SetNRelax(core::BraidCore, level::Integer, nrelax::Integer) = @ccall libbraid.braid_SetNRelax(core._braid_core::Ptr{C_NULL}, level::Cint, nrelax::Cint)::Cint
SetCRelaxWt(core::BraidCore, level::Integer, Cwt::Real) = @ccall libbraid.braid_SetCRelaxWt(core._braid_core::Ptr{C_NULL}, level::Cint, Cwt::Cdouble)::Cint
SetCFactor(core::BraidCore, level::Integer, cfactor::Integer) = @ccall libbraid.braid_SetCFactor(core._braid_core::Ptr{C_NULL}, level::Cint, cfactor::Cint)::Cint
SetMaxIter(core::BraidCore, max_iter::Integer) = @ccall libbraid.braid_SetMaxIter(core._braid_core::Ptr{C_NULL}, max_iter::Cint)::Cint
SetFMG(core::BraidCore) = @ccall libbraid.braid_SetFMG(core._braid_core::Ptr{C_NULL})::Cint
SetNFMG(core::BraidCore, k::Integer) = @ccall libbraid.braid_SetNFMG(core._braid_core::Ptr{C_NULL}, k::Cint)::Cint
SetNFMGVcyc(core::BraidCore, nfmg_Vcyc::Integer) = @ccall libbraid.braid_SetNFMGVcyc(core._braid_core::Ptr{C_NULL}, nfmg_Vcyc::Cint)::Cint
SetStorage(core::BraidCore, storage::Bool) = @ccall libbraid.braid_SetStorage(core._braid_core::Ptr{C_NULL}, storage::Cint)::Cint
SetTemporalNorm(core::BraidCore, tnorm::Bool) = @ccall libbraid.braid_SetTemporalNorm(core._braid_core::Ptr{C_NULL}, tnorm::Cint)::Cint
SetPeriodic(core::BraidCore, periodic::Bool) = @ccall libbraid.braid_SetPeriodic(core._braid_core::Ptr{C_NULL}, periodic::Cint)::Cint
SetPrintLevel(core::BraidCore, print_level::Integer) = @ccall libbraid.braid_SetPrintLevel(core._braid_core::Ptr{C_NULL}, print_level::Cint)::Cint
SetFileIOLevel(core::BraidCore, io_level::Integer) = @ccall libbraid.braid_SetFileIOLevel(core._braid_core::Ptr{C_NULL}, io_level::Cint)::Cint
SetPrintFile(core::BraidCore, filename::String) = @ccall libbraid.braid_SetPrintFile(core._braid_core::Ptr{C_NULL}, filename::Cstring)::Cint
SetDefaultPrintFile(core::BraidCore) = @ccall libbraid.braid_SetDefaultPrintFile(core._braid_core::Ptr{C_NULL})::Cint
SetAccessLevel(core::BraidCore, access_level::Integer) = @ccall libbraid.braid_SetAccessLevel(core._braid_core::Ptr{C_NULL}, access_level::Cint)::Cint
SetFinalFCRelax(core::BraidCore) = @ccall libbraid.braid_SetFinalFCRelax(core._braid_core::Ptr{C_NULL})::Cint
function GetNumIter(core::BraidCore)
    niter = Ref(0)
    @ccall libbraid.braid_GetNumIter(core._braid_core::Ptr{C_NULL}, niter::Ref{Cint})::Cint
    return niter[] # dereference
end

function GetRNorms(core::BraidCore)
    nrequest = Ref(GetNumIter(core))
    rnorms = zeros(nrequest)
    @ccall libbraid.braid_GetRNorms(core._braid_core::Ptr{C_NULL}, nrequest::Ref{Cint}, rnorms::Ref{Cdouble})::Cint
    return rnorms
end

function GetNLevels(core::BraidCore)
    nlevels = Ref(0)
    @ccall libbraid.braid_GetNLevels(core._braid_core::Ptr{C_NULL}, nlevels::Ref{Cint})::Cint
    return nlevels[]
end

SetSeqSoln(core::BraidCore, seq_soln::Bool) = @ccall libbraid.braid_GetNLevels(core._braid_core::Ptr{C_NULL}, seq_soln::Cint)::Cint

function SetRichardsonEstimation(core::BraidCore, est_error::Bool, richardson::Bool, local_order::Integer)
    @ccall libbraid.braid_SetRichardsonEstimation(core._braid_core::Ptr{C_NULL}, est_error::Cint, richardson::Cint, local_order::Cint)::Cint
end

function SetDeltaCorrection(core::BraidCore, rank::Integer, basis_init::Function, inner_prod::Function)
    core._braid_app.basis_init = basis_init
    core._braid_app.inner_prod = inner_prod
    @ccall libbraid.braid_SetDeltaCorrection(core._braid_core::Ptr{C_NULL}, rank::Cint, _c_init_basis::Ptr{Cvoid}, _c_inner_prod::Ptr{Cvoid})::Cint
end

function SetDeferDelta(core::BraidCore, level::Integer, iter::Integer)
    @ccall libbraid.braid_SetDeferDelta(core._braid_core::Ptr{C_NULL}, level::Cint, iter::Cint)::Cint
end

function SetLyapunovEstimation(core::BraidCore, relax::Bool, cglv::Bool, exponents::Bool)
    @ccall libbraid.braid_SetLyapunovEstimation(core._braid_core::Ptr{C_NULL}, relax::Cint, cglv::Cint, exponents::Cint)::Cint
end

# Still missing spatial coarsening, sync, residual, and adjoint

# braid_test
function testInitAccess(app::BraidApp, t::Real)
    @ccall libbraid.braid_TestInitAccess(
        app::Ref{BraidApp},
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

function testClone(app::BraidApp, t::Real)
    @ccall libbraid.braid_TestClone(
        app::Ref{BraidApp},
        app.comm::MPI.MPI_Comm,
        C_stdout::Base.Libc.FILE,
        t::Cdouble,
        _c_init::Ptr{Cvoid},
        _c_access::Ptr{Cvoid},
        _c_free::Ptr{Cvoid},
        _c_clone::Ptr{Cvoid}
    )::Cint
end

function testSum(app::BraidApp, t::Real)
    @ccall libbraid.braid_TestSum(
        app::Ref{BraidApp},
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

function testSpatialNorm(app::BraidApp, t::Real)
    @ccall libbraid.braid_TestSpatialNorm(
        app::Ref{BraidApp},
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

function testBuf(app::BraidApp, t::Real)
    @ccall libbraid.braid_TestBuf(
        app::Ref{BraidApp},
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

function testDelta(app::BraidApp, t::Real, dt::Real, rank::Integer)
    app.basis_init::Function
    app.inner_prod::Function
    # if app.basis_init === nothing || app.inner_prod === nothing
    #     error("User-defined functions inner_prod(u, v) and basis_init(t, index) must not be nothing.")
    # end

    @ccall libbraid.braid_TestDelta(
        app::Ref{BraidApp},
        app.comm::MPI.MPI_Comm,
        C_stdout::Base.Libc.FILE,
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
        _c_step::Ptr{Cvoid}
    )::Cint
    print('\n')
    println("Check output for objects not freed:")
    println(app.ref_ids)
end 

# END MODULE XBraid
end
