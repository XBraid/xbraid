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
using Serialization: serialize, deserialize, Serializer # used to pack arbitrary julia objects into buffers
using LinearAlgebra: norm2
using MPI
using MethodAnalysis

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

c_stdout = Libc.FILE(Libc.RawFD(1), "w")  # corresponds to C standard output
function malloc_null_double_ptr(T::Type)
	pp = Base.Libc.malloc(sizeof(Ptr{Cvoid}))
	pp = reinterpret(Ptr{Ptr{T}}, pp)
	# unsafe_store!(pp, C_NULL)
	return pp
end

# This can contain anything (TODO: do I even need this?)
mutable struct BraidVector{T}
	user_vector::T
	VecType::Type
end
BraidVector(u::T) where T = BraidVector(u, T)

OptionalFunction = Union{Function, Nothing}
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

	basis_init::OptionalFunction
	inner_prod::OptionalFunction

	# dictionary to store globally scoped references to allocated braid_vector objects
	ref_ids::IdDict{UInt64, BraidVector}   
	bufsize::Integer  # expected serialized size of user's vector
	bufsize_lyap::Integer
	user_AppType::Type
	user_VecType::Type
	user_BasType::Type
end

# some default functions assuming the user's vector is a julia array
function default_sum!(app, α::Real, x::AbstractArray, β::Real, y::AbstractArray)
	y .= α .* x .+ β .* y
	return
end

function default_norm(app, u::AbstractArray)
	return norm2(u)
end

function default_inner_prod(app, u::AbstractArray, v::AbstractArray)
	return dot(u, v)
end

# default constructor
function BraidApp(app, comm::MPI.Comm, comm_t::MPI.Comm, step::Function, init::Function, sum::Function, norm::Function, access::Function, basis_init::Function, inner_prod::Function)
	BraidApp(app, comm, comm_t, step, init, sum, norm, access, basis_init, inner_prod, IdDict(), 0, 0, Nothing, Nothing, Nothing)
end

function BraidApp(app, comm::MPI.Comm, comm_t::MPI.Comm, step::Function, init::Function, sum::Function, norm::Function, access::Function)
	BraidApp(app, comm, comm_t, step, init, sum, norm, access, nothing, nothing, IdDict(), 0, 0, Nothing, Nothing, Nothing)
end

function BraidApp(app, comm::MPI.Comm, step::Function, init::Function, sum::Function, norm::Function, access::Function)
	BraidApp(app, comm, comm, step, init, sum, norm, access)
end

function BraidApp(app, comm::MPI.Comm, step, init, access)
	BraidApp(app, comm, comm, step, init, default_sum!, default_norm, access)
end

function BraidApp(app, comm::MPI.Comm, step, init, access, basis_init)
	BraidApp(app, comm, comm, step, init, default_sum!, default_norm, access, basis_init, default_inner_prod)
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
	tstart::Real
	tstop::Real
	ntime::Integer
	function BraidCore(_braid_core, _braid_app, tstart, tstop, ntime)
		x = new(_braid_core, _braid_app, tstart, tstop, ntime)
		finalizer(x) do core
			# @async println("Destroying BraidCore $(core._braid_core)")
			@ccall libbraid.braid_Destroy(core._braid_core::Ptr{Cvoid})::Cint
		end
	end
end

"""
Used to add/remove a reference to the vector u to the IdDict stored in the app.
this keeps all newly allocated braid_vectors in the global scope, preventing them from being garbage collected. 
"""
function _register_vector(app::BraidApp, u::BraidVector)
	app.ref_ids[objectid(u)] = u
end

function _deregister_vector(app::BraidApp, u::BraidVector)
	id = objectid(u)::UInt64
	pop!(app.ref_ids, id)
end

function postInitPrecompile(app::BraidApp)
	#= 
	# This is a hack to make sure every processor knows how big the user's vector will be.
	# We can also take this time to precompile the user's step function so it doesn't happen
	# in serial on the coarse grid.
	=#
	GC.enable(false) # disable garbage collection
	app_ptr = pointer_from_objref(app)::Ptr{Cvoid}
	pp = malloc_null_double_ptr(Cvoid)

	_jl_init!(app_ptr, 0., pp)
	u_ptr = unsafe_load(pp)
	u = unsafe_pointer_to_objref(u_ptr)::BraidVector
	VecType = typeof(u.user_vector)
	AppType = typeof(app.user_app)
	
	!precompile(app.step, 	    (AppType, Ptr{Cvoid}, VecType, VecType, Float64, Float64)) && println("failed to precompile step")
	!precompile(app.spatialnorm, (AppType, VecType)) && println("failed to precompile norm")
	!precompile(app.sum,         (AppType, Float64, VecType, Float64, VecType)) && println("failed to precompile sum")
	if app.access !== nothing 
		!precompile(app.access,  (AppType, Ptr{Cvoid}, VecType)) && println("failed to precompile access")
	end

	# some julia functions that can be precompiled
	precompile(deepcopy, (BraidVector{VecType},))
	precompile(unsafe_store!, (Ptr{Float64}, Float64,))
	precompile(Tuple{typeof(Base.getproperty), BraidVector{VecType}, Symbol})
	precompile(Tuple{typeof(Base.unsafe_store!), Ptr{Int32}, Int64})
	precompile(Tuple{typeof(deserialize), Serializer{Base.GenericIOBuffer{Array{UInt8, 1}}}, DataType})
	precompile(Tuple{typeof(Base.unsafe_wrap), Type{Array{UInt8, 1}}, Ptr{UInt8}, Int64})
	precompile(Tuple{Type{NamedTuple{(:read, :write, :maxsize), T} where T<:Tuple}, Tuple{Bool, Bool, Int64}})

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

	_jl_init_basis!(app_ptr, 0., Int32(0), pp)
	ψ_ptr = unsafe_load(pp)
	ψ = unsafe_pointer_to_objref(ψ_ptr)
	BasType = typeof(ψ.user_vector)
	println("precompiling user Delta functions")
	!precompile(app.step,       (AppType, Ptr{Cvoid}, VecType, VecType, Float64, Float64, Vector{BasType})) && println("failed to compile step_du")
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
function Init(comm_world::MPI.Comm, comm_t::MPI.Comm,
	tstart::Real, tstop::Real, ntime::Integer,
	step::Function, init::Function, sum::Function, spatialnorm::Function, access::OptionalFunction;
	app = nothing
)::BraidCore
	_app = BraidApp(app, comm_world, comm_t, step, init, sum, spatialnorm, access)
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

	return BraidCore(_core, _app, tstart, tstop, ntime)
end

function Init(comm_world::MPI.Comm, tstart::Real, tstop::Real, ntime::Integer, step::Function, init::Function, access::OptionalFunction; app=nothing)
	Init(comm_world, comm_world, tstart, tstop, ntime, step, init, default_sum!, default_norm, access; app=app)
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
	fdt = (core.tstop - core.tstart) / (core.ntime+1)
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
		# println("Calling Warmup")
		Warmup(core)
	else
		_app = core._braid_app
		# println("Calling Precompile")
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


function PrintStats(core::BraidCore)
	GC.@preserve core begin
		@ccall libbraid.braid_PrintStats(core._braid_core::Ptr{Cvoid})::Cint
	end
end

function SetTimerFile(core::BraidCore, filestem::String)
	len = @ccall strlen(filestem::Cstring)::Cint
	GC.@preserve core begin
		@ccall libbraid.braid_SetTimerFile(core._braid_core::Ptr{Cvoid}, len::Cint, filestem::Cstring)::Cint
	end
end

function PrintTimers(core::BraidCore)
	GC.@preserve core begin
		@ccall libbraid.braid_PrintTimers(core._braid_core::Ptr{Cvoid})::Cint
	end
end

function ResetTimer(core::BraidCore)
	GC.@preserve core begin
		@ccall libbraid.braid_ResetTimer(core._braid_core::Ptr{Cvoid})::Cint
	end
end

function SetTimings(core::BraidCore, timing_level::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetTimings(core._braid_core::Ptr{Cvoid}, timing_level::Cint)::Cint
	end
end

function WriteConvHistory(core::BraidCore, filename::String)
	GC.@preserve core begin
		@ccall libbraid.braid_WriteConvHistory(core._braid_core::Ptr{Cvoid}, filename::Cstring)::Cint
	end
end

function SetMaxLevels(core::BraidCore, max_levels::Integer)
	@assert max_levels > 0 "Max. levels must be an integer greater than 0"
	GC.@preserve core begin
		@ccall libbraid.braid_SetMaxLevels(core._braid_core::Ptr{Cvoid}, max_levels::Cint)::Cint
	end
end

function SetIncrMaxLevels(core::BraidCore)
	GC.@preserve core begin
		@ccall libbraid.braid_SetIncrMaxLevels(core._braid_core::Ptr{Cvoid})::Cint
	end
end

function SetSkip(core::BraidCore, skip::Bool)
	GC.@preserve core begin
		@ccall libbraid.braid_SetSkip(core._braid_core::Ptr{Cvoid}, skip::Cint)::Cint
	end
end


function SetRefine(core::BraidCore, refine::Bool)
	GC.@preserve core begin
		@ccall libbraid.braid_SetSkip(core._braid_core::Ptr{Cvoid}, refine::Cint)::Cint
	end
end

function SetMaxRefinements(core::BraidCore, max_refinements::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetMaxRefinements(core._braid_core::Ptr{Cvoid}, max_refinements::Cint)::Cint
	end
end


function SetTPointsCutoff(core::BraidCore, tpoints_cutoff::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetTPointsCutoff(core._braid_core::Ptr{Cvoid}, tpoints_cutoff::Cint)::Cint
	end
end


function SetMinCoarse(core::BraidCore, min_coarse::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetMinCoarse(core._braid_core::Ptr{Cvoid}, min_coarse::Cint)::Cint
	end
end

function SetRelaxOnlyCG(core::BraidCore, relax_only_cg::Bool)
	GC.@preserve core begin
		@ccall libbraid.braid_SetRelaxOnlyCG(core._braid_core::Ptr{Cvoid}, relax_only_cg::Cint)::Cint
	end
end

function SetAbsTol(core::BraidCore, atol::Real)
	GC.@preserve core begin
		@ccall libbraid.braid_SetAbsTol(core._braid_core::Ptr{Cvoid}, atol::Cdouble)::Cint
	end
end

function SetRelTol(core::BraidCore, rtol::Real)
	GC.@preserve core begin
		@ccall libbraid.braid_SetRelTol(core._braid_core::Ptr{Cvoid}, rtol::Cdouble)::Cint
	end
end

function SetNRelax(core::BraidCore, level::Integer, nrelax::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetNRelax(core._braid_core::Ptr{Cvoid}, level::Cint, nrelax::Cint)::Cint
	end
end

function SetCRelaxWt(core::BraidCore, level::Integer, Cwt::Real)
	GC.@preserve core begin
		@ccall libbraid.braid_SetCRelaxWt(core._braid_core::Ptr{Cvoid}, level::Cint, Cwt::Cdouble)::Cint
	end
end

function SetCFactor(core::BraidCore, level::Integer, cfactor::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetCFactor(core._braid_core::Ptr{Cvoid}, level::Cint, cfactor::Cint)::Cint
	end
end

function SetMaxIter(core::BraidCore, max_iter::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetMaxIter(core._braid_core::Ptr{Cvoid}, max_iter::Cint)::Cint
	end
end

function SetFMG(core::BraidCore)
	GC.@preserve core begin
		@ccall libbraid.braid_SetFMG(core._braid_core::Ptr{Cvoid})::Cint
	end
end

function SetNFMG(core::BraidCore, k::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetNFMG(core._braid_core::Ptr{Cvoid}, k::Cint)::Cint
	end
end

function SetNFMGVcyc(core::BraidCore, nfmg_Vcyc::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetNFMGVcyc(core._braid_core::Ptr{Cvoid}, nfmg_Vcyc::Cint)::Cint
	end
end

function SetStorage(core::BraidCore, storage::Bool)
	GC.@preserve core begin
		@ccall libbraid.braid_SetStorage(core._braid_core::Ptr{Cvoid}, storage::Cint)::Cint
	end
end

function SetTemporalNorm(core::BraidCore, tnorm::Bool)
	GC.@preserve core begin
		@ccall libbraid.braid_SetTemporalNorm(core._braid_core::Ptr{Cvoid}, tnorm::Cint)::Cint
	end
end

function SetPeriodic(core::BraidCore, periodic::Bool)
	GC.@preserve core begin
		@ccall libbraid.braid_SetPeriodic(core._braid_core::Ptr{Cvoid}, periodic::Cint)::Cint
	end
end

function SetPrintLevel(core::BraidCore, print_level::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetPrintLevel(core._braid_core::Ptr{Cvoid}, print_level::Cint)::Cint
	end
end

function SetFileIOLevel(core::BraidCore, io_level::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetFileIOLevel(core._braid_core::Ptr{Cvoid}, io_level::Cint)::Cint
	end
end

function SetPrintFile(core::BraidCore, filename::String)
	GC.@preserve core begin
		@ccall libbraid.braid_SetPrintFile(core._braid_core::Ptr{Cvoid}, filename::Cstring)::Cint
	end
end

function SetDefaultPrintFile(core::BraidCore)
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
function SetAccessLevel(core::BraidCore, access_level::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetAccessLevel(core._braid_core::Ptr{Cvoid}, access_level::Cint)::Cint
	end
end

function SetFinalFCRelax(core::BraidCore)
	GC.@preserve core begin
		@ccall libbraid.braid_SetFinalFCRelax(core._braid_core::Ptr{Cvoid})::Cint
	end
end

function GetNumIter(core::BraidCore)
	niter = Ref{Cint}(0)
	GC.@preserve core begin
		@ccall libbraid.braid_GetNumIter(core._braid_core::Ptr{Cvoid}, niter::Ref{Cint})::Cint
	end
	return niter[]
end

function GetRNorms(core::BraidCore)
	nrequest = Ref{Cint}()
	nrequest[] = GetNumIter(core)
	rnorms = zeros(nrequest[])
	GC.@preserve core begin
		@ccall libbraid.braid_GetRNorms(core._braid_core::Ptr{Cvoid}, nrequest::Ref{Cint}, rnorms::Ref{Cdouble})::Cint
	end
	nrequest[] == 0 && return Float64[]
	return rnorms[1:nrequest[]]
end

function GetNLevels(core::BraidCore)
	nlevels = Ref(0)
	GC.@preserve core begin
		@ccall libbraid.braid_GetNLevels(core._braid_core::Ptr{Cvoid}, nlevels::Ref{Cint})::Cint
	end
	return nlevels[]
end

function SetSeqSoln(core::BraidCore, seq_soln::Bool)
	GC.@preserve core begin
		@ccall libbraid.braid_GetNLevels(core._braid_core::Ptr{Cvoid}, seq_soln::Cint)::Cint
	end
end

function SetRichardsonEstimation(core::BraidCore, est_error::Bool, richardson::Bool, local_order::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetRichardsonEstimation(core._braid_core::Ptr{Cvoid}, est_error::Cint, richardson::Cint, local_order::Cint)::Cint
	end
end

function SetDeltaCorrection(core::BraidCore, rank::Integer, basis_init::Function, inner_prod::Function)
	core._braid_app.basis_init = basis_init
	core._braid_app.inner_prod = inner_prod
	# deltaPrecompile(core._braid_app)

	GC.@preserve core begin
		@ccall libbraid.braid_SetDeltaCorrection(core._braid_core::Ptr{Cvoid}, rank::Cint, _c_init_basis::Ptr{Cvoid}, _c_inner_prod::Ptr{Cvoid})::Cint
	end
end

function SetDeferDelta(core::BraidCore, level::Integer, iter::Integer)
	GC.@preserve core begin
		@ccall libbraid.braid_SetDeferDelta(core._braid_core::Ptr{Cvoid}, level::Cint, iter::Cint)::Cint
	end
end

function SetLyapunovEstimation(core::BraidCore; relax::Bool = false, cglv::Bool = true, exponents::Bool = false)
	GC.@preserve core begin
		@ccall libbraid.braid_SetLyapunovEstimation(core._braid_core::Ptr{Cvoid}, relax::Cint, cglv::Cint, exponents::Cint)::Cint
	end
end

# Still missing spatial coarsening, sync, residual, and adjoint

# braid_test
function testInitAccess(app::BraidApp, t::Real, outputFile::Libc.FILE)
	@ccall libbraid.braid_TestInitAccess(
		app::Ref{BraidApp},
		app.comm::MPI.MPI_Comm,
		outputFile::Libc.FILE,
		t::Cdouble,
		_c_init::Ptr{Cvoid},
		_c_access::Ptr{Cvoid},
		_c_free::Ptr{Cvoid},
	)::Cint

	# println("Serialized size of user vector: $(app.bufsize)")
	# println("Check output for objects not properly freed:")
	# println(app.ref_ids)
end
testInitAccess(app::BraidApp, t::Real, outputFile::IO) = testInitAccess(app, t, Libc.FILE(outputFile))
testInitAccess(app::BraidApp, t::Real) = testInitAccess(app, t, c_stdout)

function testClone(app::BraidApp, t::Real, outputFile::Libc.FILE)
	@ccall libbraid.braid_TestClone(
		app::Ref{BraidApp},
		app.comm::MPI.MPI_Comm,
		outputFile::Libc.FILE,
		t::Cdouble,
		_c_init::Ptr{Cvoid},
		_c_access::Ptr{Cvoid},
		_c_free::Ptr{Cvoid},
		_c_clone::Ptr{Cvoid},
	)::Cint
end
testClone(app::BraidApp, t::Real, outputFile::IO) = testClone(app, t, Libc.FILE(outputFile))
testClone(app::BraidApp, t::Real) = testClone(app, t, c_stdout)

function testSum(app::BraidApp, t::Real, outputFile::Libc.FILE)
	@ccall libbraid.braid_TestSum(
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
end
testSum(app::BraidApp, t::Real, outputFile::IO) = testSum(app, t, Libc.FILE(outputFile))
testSum(app::BraidApp, t::Real) = testSum(app, t, c_stdout)

function testSpatialNorm(app::BraidApp, t::Real, outputFile::Libc.FILE)
	@ccall libbraid.braid_TestSpatialNorm(
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
end
testSpatialNorm(app::BraidApp, t::Real, outputFile::IO) = testSpatialNorm(app, t, Libc.FILE(outputFile))
testSpatialNorm(app::BraidApp, t::Real) = testSpatialNorm(app, t, c_stdout)

function testBuf(app::BraidApp, t::Real, outputFile::Libc.FILE)
	@ccall libbraid.braid_TestBuf(
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
	# print('\n')
	# println("Check output for objects not freed:")
	# println(app.ref_ids)
end
testBuf(app::BraidApp, t::Real, outputFile::IO) = testBuf(app, t, Libc.FILE(outputFile))
testBuf(app::BraidApp, t::Real) = testBuf(app, t, c_stdout)

function testDelta(app::BraidApp, t::Real, dt::Real, rank::Integer, outputFile::Libc.FILE)
	app.basis_init::Function
	app.inner_prod::Function
	@ccall libbraid.braid_TestDelta(
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
	print('\n')
	# println("Check output for objects not freed:")
	# println(app.ref_ids)
end
testDelta(app::BraidApp, t::Real, dt::Real, rank::Integer, outputFile::IO) = testDelta(app, t, dt, rank, Libc.FILE(outputFile))
testDelta(app::BraidApp, t::Real, dt::Real, rank::Integer) = testDelta(app, t, dt, rank, c_stdout)

# END MODULE XBraid
end
