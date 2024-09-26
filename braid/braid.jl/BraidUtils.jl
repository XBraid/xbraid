module BraidUtils
export libbraid, c_stdout, BlackHoleBuffer, malloc_null_double_ptr, stacktrace_warn

libbraid = joinpath(dirname(@__DIR__), "libbraid.so")
c_stdout = Libc.FILE(Libc.RawFD(1), "w")  # corresponds to C standard output


"""
Allocates a valid pointer to a pointer of type T
"""
function malloc_null_double_ptr(T::Type)
	pp = Base.Libc.malloc(sizeof(Ptr{Cvoid}))
	pp = reinterpret(Ptr{Ptr{T}}, pp)
	return pp
end

"""
Displays a caught error message, including stacktrace, without throwing
"""
function stacktrace_warn(msg::String, err)
    err_msg = sprint(showerror, err)
    trace = sprint((io,v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
    @warn "$(msg):\n$(err_msg)\n$(trace)"
end

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
    return 1
end
function Base.write(to::BlackHoleBuffer, x::Array{T}) where T
    to.ptr += sizeof(x)
    return sizeof(x)
end

end # module BraidUtils

# helper functions
isCPoint(i::Integer, cfactor::Integer)::Bool = ((i-1) % cfactor == 0)
isFPoint(i::Integer, cfactor::Integer)::Bool = !isCPoint(i, cfactor)
mapFineToCoarse(i::Integer, cfactor::Integer)::Integer = (i-1) ÷ cfactor + 1
mapCoarseToFine(i::Integer, cfactor::Integer)::Integer = (i-1) * cfactor + 1