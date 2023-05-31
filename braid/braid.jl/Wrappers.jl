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

module Wrappers

using LinearAlgebra: norm2, dot
using Serialization: serialize, deserialize
using ..XBraid.BraidUtils
using ..XBraid: Status, BraidApp, BraidVector

# useful helpers

"""
Macro which wraps the function call in a try-catch block which prints the stacktrace without throwing,
since Julia functions should not throw errors when called from braid
"""
macro braidWrapper(expr)
    # this error is at compile time, and therefore fine
    expr.head == :function || error("Not a function expression.")
    funcname = expr.args[1].args[1].args[1]
    user_funcname = split(string(funcname), "_")[end]
    :(
        function $(esc(funcname))(args...; kwargs...)::Cint
            f = $expr
            try
                result = f(args...; kwargs...)
                return result
            catch err
                stacktrace_warn("Error in user function " * $(esc(user_funcname)), err)
                return 1
            end
        end
    )
end

"""
Used to add a reference to the vector u to the IdDict stored in the app.
this keeps all newly allocated braid_vectors in the global scope, preventing them from being garbage collected. 
"""
function _register_vector(app::BraidApp, u::BraidVector)
    app.ref_ids[objectid(u)] = u
end

"""
Used to remove a reference to the vector u to the IdDict stored in the app.
"""
function _deregister_vector(app::BraidApp, u::BraidVector)
    id = objectid(u)::UInt64
    pop!(app.ref_ids, id)
end

# these are internal functions which directly interface with XBraid


# Step, Init, and Access are required

@braidWrapper function _jl_step!(_app::Ptr{Cvoid},
    _ustop::Ptr{Cvoid},
    _fstop::Ptr{Cvoid},
    _u::Ptr{Cvoid},
    _status::Ptr{Cvoid}
)::Cint
    # println("step")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    ustop = unsafe_pointer_to_objref(_ustop)::BraidVector
    status = Status.StepStatus(_status)

    tstart, tstop = Status.getTstartTstop(status)
    delta_rank = Status.getDeltaRank(status)

    # call the user's function
    # residual option
    if _fstop !== C_NULL
        fstop = unsafe_pointer_to_objref(_fstop)::BraidVector
        app.step(app.user_app, status, u.user_vector, ustop.user_vector, fstop.user_vector, tstart, tstop)
        # Delta correction
    elseif delta_rank > 0
        basis_vecs = Status.getBasisVectors(status)
        app.step(app.user_app, status, u.user_vector, ustop.user_vector, tstart, tstop, basis_vecs)
        # Default
    else
        app.step(app.user_app, status, u.user_vector, ustop.user_vector, tstart, tstop)
    end

    return 0
end
precompile(_jl_step!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
_c_step = @cfunction(_jl_step!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
export _c_step, _jl_step!

@braidWrapper function _jl_init!(_app::Ptr{Cvoid}, t::Cdouble, u_ptr::Ptr{Ptr{Cvoid}})::Cint
    # println("init")
    app = unsafe_pointer_to_objref(_app)::BraidApp

    # initialize u and register a reference with IdDict
    u = BraidVector(app.init(app.user_app, t))
    _register_vector(app, u)

    unsafe_store!(u_ptr, pointer_from_objref(u))

    # store max size of all initialized vectors
    if app.bufsize == 0
        # This serializes u but doesn't actually
        # store the data
        buffer = BlackHoleBuffer()
        serialize(buffer, u)
        app.bufsize = buffer.ptr + sizeof(Int)
    end

    if app.user_VecType == Nothing
        app.user_VecType = u.VecType
    end

    return 0
end
_c_init = @cfunction(_jl_init!, Cint, (Ptr{Cvoid}, Cdouble, Ptr{Ptr{Cvoid}}))
export _c_init, _jl_init!

@braidWrapper function _jl_init_basis!(_app::Ptr{Cvoid}, t::Cdouble, index::Cint, u_ptr::Ptr{Ptr{Cvoid}})::Cint
    # println("init_basis")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    # julia uses 1-based indexing
    u = BraidVector(app.basis_init(app.user_app, t, index + 1))

    _register_vector(app, u)
    unsafe_store!(u_ptr, pointer_from_objref(u))

    # store max size of all initialized vectors
    if app.bufsize == 0
        buffer = BlackHoleBuffer()
        serialize(buffer, u)
        app.bufsize_lyap = buffer.ptr + sizeof(Int)
    end

    # store type of initialized vector
    if app.user_BasType == Nothing
        app.user_BasType = u.VecType
    end

    return 0
end
_c_init_basis = @cfunction(_jl_init_basis!, Cint, (Ptr{Cvoid}, Cdouble, Cint, Ptr{Ptr{Cvoid}}))
export _c_init_basis, _jl_init_basis!

@braidWrapper function _jl_clone!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, v_ptr::Ptr{Ptr{Cvoid}})::Cint
    # println("clone")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    # initialize v, and copy u into v
    v = deepcopy(u)

    # then register v with IdDict and store in v_ptr
    _register_vector(app, v)
    unsafe_store!(v_ptr, pointer_from_objref(v))

    return 0
end
precompile(_jl_clone!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}))
_c_clone = @cfunction(_jl_clone!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}))
export _c_clone, _jl_clone!

@braidWrapper function _jl_free!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid})::Cint
    # println("free")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    # removing the global reference to u will cause it to be garbage collected
    _deregister_vector(app, u)
    return 0
end
precompile(_jl_free!, (Ptr{Cvoid}, Ptr{Cvoid}))
_c_free = @cfunction(_jl_free!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
export _c_free, _jl_free!

@braidWrapper function _jl_sum!(_app::Ptr{Cvoid},
                                alpha::Cdouble, _x::Ptr{Cvoid},
                                beta::Cdouble, _y::Ptr{Cvoid})::Cint
    app = unsafe_pointer_to_objref(_app)::BraidApp
    x = unsafe_pointer_to_objref(_x)::BraidVector
    y = unsafe_pointer_to_objref(_y)::BraidVector
    app.sum(app.user_app, alpha, x.user_vector, beta, y.user_vector)
    return 0
end
precompile(_jl_sum!, (Ptr{Cvoid}, Cdouble, Ptr{Cvoid}, Cdouble, Ptr{Cvoid}))
_c_sum = @cfunction(_jl_sum!, Cint, (Ptr{Cvoid}, Cdouble, Ptr{Cvoid}, Cdouble, Ptr{Cvoid}))
export _c_sum, _jl_sum!

@braidWrapper function _jl_norm!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, norm_ptr::Ptr{Cdouble})::Cint
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    norm = app.spatialnorm(app.user_app, u.user_vector)
    unsafe_store!(norm_ptr, norm)

    return 0
end
precompile(_jl_norm!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))
_c_norm = @cfunction(_jl_norm!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))
export _c_norm, _jl_norm!

@braidWrapper function _jl_inner_prod!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, _v::Ptr{Cvoid}, norm_ptr::Ptr{Cdouble})::Cint
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    v = unsafe_pointer_to_objref(_v)::BraidVector
    prod = app.inner_prod(app.user_app, u.user_vector, v.user_vector)
    unsafe_store!(norm_ptr, prod)

    return 0
end
precompile(_jl_inner_prod!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))
_c_inner_prod = @cfunction(_jl_inner_prod!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))
export _c_inner_prod, _jl_inner_prod!

@braidWrapper function _jl_access!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, _status::Ptr{Cvoid})::Cint
    app = unsafe_pointer_to_objref(_app)::BraidApp
    status = Status.AccessStatus(_status)

    if !isnothing(app.access)
        u = unsafe_pointer_to_objref(_u)::BraidVector
        app.access(app.user_app, status, u.user_vector)
    end

    return 0
end
precompile(_jl_access!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
_c_access = @cfunction(_jl_access!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
export _c_access, _jl_access!

# buffer functions

@braidWrapper function _jl_bufsize!(_app::Ptr{Cvoid}, size_ptr::Ptr{Cint}, status::Ptr{Cvoid})::Cint
    app = unsafe_pointer_to_objref(_app)
    unsafe_store!(size_ptr, app.bufsize)
    return 0
end
_c_bufsize = @cfunction(_jl_bufsize!, Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cvoid}))
export _c_bufsize, _jl_bufsize!

@braidWrapper function _jl_bufpack!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, _buffer::Ptr{Cvoid}, status::Ptr{Cvoid})::Cint
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    buff_arr = unsafe_wrap(Vector{UInt8}, Base.unsafe_convert(Ptr{UInt8}, _buffer), app.bufsize)
    buffer = IOBuffer(buff_arr, read=true, write=true, maxsize=app.bufsize)

    # store u in buffer
    seek(buffer, sizeof(Int))
    serialize(buffer, u)

    # also store the total size of the buffer
    seek(buffer, 0)
    write(buffer, buffer.size)

    # tell braid the written size
    @ccall libbraid.braid_BufferStatusSetSize(status::Ptr{Cvoid}, buffer.size::Cdouble)::Cint

    return 0
end
precompile(_jl_bufpack!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
_c_bufpack = @cfunction(_jl_bufpack!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
export _c_bufpack, _jl_bufpack!

@braidWrapper function _jl_bufunpack!(_app::Ptr{Cvoid}, _buffer::Ptr{Cvoid}, u_ptr::Ptr{Ptr{Cvoid}}, status::Ptr{Cvoid})::Cint
    @assert _buffer !== C_NULL "tried to unpack null buffer"
    app = unsafe_pointer_to_objref(_app)::BraidApp
    # get size of buffer we are unpacking:
    header = reinterpret(Ptr{Int}, _buffer)
    bufsize = unsafe_load(header)::Int
    buff_arr = unsafe_wrap(Vector{UInt8}, Base.unsafe_convert(Ptr{UInt8}, _buffer), bufsize)
    buffer = IOBuffer(buff_arr, read=true, write=true, maxsize=bufsize)
    seek(buffer, sizeof(Int))

    # unpack the buffer into a new julia object, then register with IdDict
    u = deserialize(buffer)::BraidVector{app.user_VecType}
    _register_vector(app, u)

    # store u in provided pointer
    unsafe_store!(u_ptr, pointer_from_objref(u))
    return 0
end
precompile(_jl_bufunpack!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
_c_bufunpack = @cfunction(_jl_bufunpack!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
export _c_bufunpack, _jl_bufunpack!

# optional functions

@braidWrapper function _jl_sync(_app::Ptr{Cvoid}, status::Ptr{Cvoid})::Cint
    app = unsafe_pointer_to_objref(_app)::BraidApp
    if app.sync !== nothing
        app.sync(app.user_app, Status.SyncStatus(status))
    end
    return 0
end
precompile(_jl_sync, (Ptr{Cvoid}, Ptr{Cvoid}))
_c_sync = @cfunction(_jl_sync, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
export _c_sync, _jl_sync

# default functions assuming the user's vector is a julia array
export default_sum!, default_norm, default_inner_prod

function default_sum!(app, α::Real, x::AbstractArray, β::Real, y::AbstractArray)
    y .= α .* x .+ β .* y
end
default_norm(app, u::AbstractArray) = norm2(u)
default_inner_prod(app, u::AbstractArray, v::AbstractArray) = dot(u, v)


end # module Wrappers