# TODO: figure out status wrapper
# @ccall requires these type annotations to be known at compile time, so macros won't work here
function status_GetT(status)
    t = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetT(status::Ptr{Cvoid}, t::Ref{Cdouble})::Cint
    return t[]
end

function status_GetTstartTstop(status::Ptr{Cvoid})
    tstart, tstop = Ref{Cdouble}(), Ref{Cdouble}()
    @ccall libbraid.braid_StepStatusGetTstartTstop(status::Ptr{Cvoid}, tstart::Ref{Cdouble}, tstop::Ref{Cdouble})::Cint
    return tstart[], tstop[]
end

function status_GetTIndex(status)
    ti = Ref{Cint}()
    @ccall libbraid.braid_StatusGetTIndex(status::Ptr{Cvoid}, ti::Ref{Cint})::Cint
    return ti[]
end

function status_GetLevel(status)
    level = Ref{Cint}()
    @ccall libbraid.braid_StatusGetLevel(status::Ptr{Cvoid}, level::Ref{Cint})::Cint
    return level[]
end

function status_GetWrapperTest(status)
    wtest = Ref{Cint}()
    @ccall libbraid.braid_StatusGetWrapperTest(status::Ptr{Cvoid}, wtest::Ref{Cint})::Cint
    return (wtest[] > 0)
end

function status_GetDeltaRank(status)
    rank = Ref{Cint}()
    @ccall libbraid.braid_StatusGetDeltaRank(status::Ptr{Cvoid}, rank::Ref{Cint})::Cint
    return rank[]
end

# These are special cases
function status_GetLocalLyapExponents(status)
    rank = status_GetDeltaRank(status)
    rank < 1 && return []

    exps = zeros(rank)
    num_retrieved = Ref(rank)
    @ccall  libbraid.braid_StatusGetLocalLyapExponents(status::Ptr{Cvoid}, exps::Ref{Cdouble}, num_retrieved::Ref{Cint})::Cint
    return exps
end

function status_GetBasisVectors(status)
    rank = status_GetDeltaRank(status)
    Ψ = []
    rank < 1 && return Ψ

    for i in 1:rank
        pp = malloc_null_double_ptr(Cvoid)
        GC.@preserve pp begin
            @ccall libbraid.braid_StatusGetBasisVec(status::Ptr{Cvoid}, pp::Ptr{Ptr{Cvoid}}, (i-1)::Cint)::Cint
        end
        p = unsafe_load(pp)
        if p !== C_NULL
            φ = unsafe_pointer_to_objref(p)::BraidVector
            push!(Ψ, φ.user_vector)
        end
        Base.Libc.free(pp)
    end

    return Ψ
end
