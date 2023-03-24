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

function status_GetIter(status)
    iter = Ref{Cint}()
    @ccall libbraid.braid_StatusGetIter(status::Ptr{Cvoid}, iter::Ref{Cint})::Cint
    return iter[]
end

function status_GetNLevels(status)
    nlevels = Ref{Cint}()
    @ccall libbraid.braid_StatusGetNLevels(status::Ptr{Cvoid}, nlevels::Ref{Cint})::Cint
    return nlevels[]
end

function status_GetNTPoints(status)
    ntpoints = Ref{Cint}()
    @ccall libbraid.braid_StatusGetNTPoints(status::Ptr{Cvoid}, ntpoints::Ref{Cint})::Cint
    return ntpoints[]
end

function status_GetResidual(status)
    res = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetResidual(status::Ptr{Cvoid}, res::Ref{Cdouble})::Cint
    return res[]
end

function status_GetDone(status)
    done = Ref{Cint}()
    @ccall libbraid.braid_StatusGetDone(status::Ptr{Cvoid}, done::Ref{Cint})::Cint
    return (done[] > 0)
end

function status_GetTILD(status)
    time = Ref{Cdouble}()
    iter = Ref{Cint}()
    level = Ref{Cint}()
    done = Ref{Cint}()
    @ccall libbraid.braid_StatusGetTILD(status::Ptr{Cvoid}, time::Ref{Cdouble}, iter::Ref{Cint}, level::Ref{Cint}, done::Ref{Cint})::Cint
    return time[], iter[], level[], (done[] > 0)
end

function status_GetCallingFunction(status)
    func = Ref{Cint}()
    @ccall libbraid.braid_StatusGetCallingFunction(status::Ptr{Cvoid}, func::Ref{Cint})::Cint
    return func[]
end

function status_GetTol(status)
    tol = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetTol(status::Ptr{Cvoid}, tol::Ref{Cdouble})::Cint
    return tol[]
end

function status_GetRNorms(status)
    iters = status_GetIter(status)
    nrequest = Ref{Cint}(iters)
    norms = zeros(Cdouble, iters)
    @ccall libbraid.braid_StatusGetRNorms(status::Ptr{Cvoid}, nrequest::Ref{Cint}, norms::Ref{Cdouble})::Cint
    return norms
end

function status_GetProc(status)
    proc = Ref{Cint}()
    @ccall libbraid.braid_StatusGetProc(status::Ptr{Cvoid}, proc::Ref{Cint})::Cint
    return proc[]
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
