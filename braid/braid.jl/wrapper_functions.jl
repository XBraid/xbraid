
"""
Where'd all my data go?
This extends Base to include a buffer which throws away the data written to it
(useful for measuring the serialized size of an object)
"""
mutable struct BlackHoleBuffer <: IO
    ptr::Int
end
BlackHoleBuffer() = BlackHoleBuffer(0)

function Base.read(from::BlackHoleBuffer, T::Type{UInt8})
    throw(ArgumentError("BlackHoleBuffer is not readable)"))
end
function Base.write(to::BlackHoleBuffer, x::UInt8)
    to.ptr += 1
    return sizeof(UInt8)
end

# these are internal functions which directly interface with XBraid
# TODO: get the status structures working
# Do we want to unpack all possible values it could contain into a julia struct?
# Or do we want to pass a pointer to the user functions that they can call StatusGet on?

function _jl_step!(_app::Ptr{Cvoid},
                   _ustop::Ptr{Cvoid},
                   _fstop::Ptr{Cvoid},
                   _u::Ptr{Cvoid},
                   status::Ptr{Cvoid})::Cint
    # println("step")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    ustop = unsafe_pointer_to_objref(_ustop)::BraidVector
    tstart, tstop = status_GetTstartTstop(status)
    delta_rank = status_GetDeltaRank(status)

    # call the user's function
    try
        if _fstop !== C_NULL
            fstop = unsafe_pointer_to_objref(_fstop)::BraidVector
            app.step(app.user_app, status, u.user_vector, ustop.user_vector, fstop.user_vector, tstart, tstop)
        elseif delta_rank > 0
            basis_vecs = status_GetBasisVectors(status)
            app.step(app.user_app, status, u.user_vector, ustop.user_vector, tstart, tstop, basis_vecs)
        else
            app.step(app.user_app, status, u.user_vector, ustop.user_vector, tstart, tstop)
        end
    catch err
        # we don't want julia exceptions to cause XBraid to exit in an undefined state...
        # so just print the message and move on
        print("Error in user Step: ")
        # println(err)
        throw(err)
    end
    # if _fstop !== C_NULL
    #     fstop = unsafe_pointer_to_objref(_fstop)::BraidVector
    #     app.step(app.user_app, status, u.user_vector, ustop.user_vector, fstop.user_vector, tstart[], tstop[])
    # elseif delta_rank[] > 0
    #     app.step(app.user_app, status, u.user_vector, ustop.user_vector, tstart[], tstop[], basis_vecs)
    # else
    #     app.step(app.user_app, status, u.user_vector, ustop.user_vector, tstart[], tstop[])
    # end

    return 0
end
precompile(_jl_step!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
_c_step = @cfunction(_jl_step!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))


function _jl_init!(_app::Ptr{Cvoid}, t::Cdouble, u_ptr::Ptr{Ptr{Cvoid}})::Cint
    # println("init")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    # initialize u and register a reference with IdDict
    u = nothing
    try
        u = BraidVector(app.init(app.user_app, t))
    catch err
        print("Error in user Init: ")
        # println(err)
        throw(err)
        return 1
    end
    # u = BraidVector(app.init(app.user_app, t))
    _register_vector(app, u)

    unsafe_store!(u_ptr, pointer_from_objref(u))

    # store max size of all initialized vectors
    if app.bufsize == 0
        # This serializes u but doesn't actually
        # store the data
        # buffer = BlackHoleBuffer()
        buffer = IOBuffer()
        serialize(buffer, u)
        app.bufsize = buffer.ptr + sizeof(Int)
    end

    if app.user_VecType == Nothing
        app.user_VecType = u.VecType
    end

    # TODO: figure out how to put an upper bound on the size of the serialized object
    # without serializing it first
    # u_size = Base.summarysize(Ref(u)) + 9
    # if u_size > app.bufsize
    #     app.bufsize = u_size
    # end

    return 0
end
_c_init = @cfunction(_jl_init!, Cint, (Ptr{Cvoid}, Cdouble, Ptr{Ptr{Cvoid}}))

function _jl_init_basis!(_app::Ptr{Cvoid}, t::Cdouble, index::Cint, u_ptr::Ptr{Ptr{Cvoid}})::Cint
    # println("init_basis")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    # try
    #     u = BraidVector(app.basis_init(app.user_app, t, index))
    # catch err
    #     print("Error in user InitBasis: ")
    #     println(err)
    #     return 1
    # end
    u = BraidVector(app.basis_init(app.user_app, t, index))
    _register_vector(app, u)
    unsafe_store!(u_ptr, pointer_from_objref(u))

    # store max size of all initialized vectors
    if app.bufsize == 0
        buffer = IOBuffer()
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

function _jl_clone!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, v_ptr::Ptr{Ptr{Cvoid}})::Cint
    # println("clone")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    # initialize v, and copy u into v
    v = nothing
    try
        v = deepcopy(u)
    catch
        return 1
    end

    # then register v with IdDict and store in v_ptr
    _register_vector(app, v)
    unsafe_store!(v_ptr, pointer_from_objref(v))

    return 0
end
precompile(_jl_clone!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}))
_c_clone = @cfunction(_jl_clone!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}))

function _jl_free!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid})::Cint
    # println("free")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    # removing the global reference to u will cause it to be garbage collected
    _deregister_vector(app, u)
    return 0
end
precompile(_jl_free!, (Ptr{Cvoid}, Ptr{Cvoid}))
_c_free = @cfunction(_jl_free!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))

function _jl_sum!(_app::Ptr{Cvoid},
    alpha::Cdouble, _x::Ptr{Cvoid},
    beta::Cdouble, _y::Ptr{Cvoid})::Cint
    # println("sum")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    x = unsafe_pointer_to_objref(_x)::BraidVector
    y = unsafe_pointer_to_objref(_y)::BraidVector
    try
        app.sum(app.user_app, alpha, x.user_vector, beta, y.user_vector)
    catch err
        print("Error in user Sum: ")
        # println(err)
        throw(err)
    end
    # app.sum(app.user_app, alpha, x.user_vector, beta, y.user_vector)
    return 0
end
precompile(_jl_sum!, (Ptr{Cvoid}, Cdouble, Ptr{Cvoid}, Cdouble, Ptr{Cvoid}))
_c_sum = @cfunction(_jl_sum!, Cint, (Ptr{Cvoid}, Cdouble, Ptr{Cvoid}, Cdouble, Ptr{Cvoid}))

function _jl_norm!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, norm_ptr::Ptr{Cdouble})::Cint
    # println("norm")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    norm = NaN
    try
        norm = app.spatialnorm(app.user_app, u.user_vector)
    catch err
        print("Error in user Norm: ")
        throw(err)
        # println(err)
    end
    # norm = app.spatialnorm(app.user_app, u.user_vector)
    unsafe_store!(norm_ptr, norm)

    return 0
end
precompile(_jl_norm!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))
_c_norm = @cfunction(_jl_norm!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))

function _jl_inner_prod!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, _v::Ptr{Cvoid}, norm_ptr::Ptr{Cdouble})::Cint
    # println("inner_prod")
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    v = unsafe_pointer_to_objref(_v)::BraidVector
    prod = NaN
    try
        prod = app.inner_prod(app.user_app, u.user_vector, v.user_vector)
    catch err
        print("Error in user InnerProd: ")
        # println(err)
        throw(err)
    end
    # prod = app.inner_prod(app.user_app, u.user_vector, v.user_vector)
    unsafe_store!(norm_ptr, prod)

    return 0
end
precompile(_jl_inner_prod!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))
_c_inner_prod = @cfunction(_jl_inner_prod!, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}))

function _jl_access!(_app::Ptr{Cvoid}, _u::Ptr{Cvoid}, status::Ptr{Cvoid})::Cint
    # println("access")
    app = unsafe_pointer_to_objref(_app)::BraidApp

    if !isnothing(app.access)
        u = unsafe_pointer_to_objref(_u)::BraidVector
        try
            app.access(app.user_app, status, u.user_vector)
        catch err
            print("Error in user Access: ")
            # println(err)
            throw(err)
        end
        # app.access(app.user_app, status, u.user_vector)
    end

    return 0
end
precompile(_jl_access!, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
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
    app = unsafe_pointer_to_objref(_app)::BraidApp
    u = unsafe_pointer_to_objref(_u)::BraidVector
    # buffer = wrap_buffer(_buffer, app.bufsize + sizeof(Int))
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

function _jl_bufunpack!(_app::Ptr{Cvoid}, _buffer::Ptr{Cvoid}, u_ptr::Ptr{Ptr{Cvoid}}, status::Ptr{Cvoid})::Cint
    # println("bufunpack")
    @assert _buffer !== C_NULL "tried to unpack null buffer"
    app = unsafe_pointer_to_objref(_app)::BraidApp
    # get size of buffer we are unpacking:
    header = reinterpret(Ptr{Int}, _buffer)
    bufsize = unsafe_load(header)::Int
    buff_arr = unsafe_wrap(Vector{UInt8}, Base.unsafe_convert(Ptr{UInt8}, _buffer), bufsize)
    buffer = IOBuffer(buff_arr, read=true, write=true, maxsize=bufsize)
    # buffer = wrap_buffer(_buffer, bufsize)
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
