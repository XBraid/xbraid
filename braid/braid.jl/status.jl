# should these automatically qeuery the braid status for some commonly used values?
struct StepStatus
    ptr::Ptr{Cvoid}
end

struct AccessStatus
    ptr::Ptr{Cvoid}
end

# @ccall requires these type annotations to be known at compile time, so macros won't work here
function status_GetT(status::Union{StepStatus, AccessStatus})
    t = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetT(status.ptr::Ptr{Cvoid}, t::Ref{Cdouble})::Cint
    return t[]
end

function status_GetTStop(status::StepStatus)
    tstop = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetTStop(status.ptr::Ptr{Cvoid}, tstop::Ref{Cdouble})::Cint
    return tstop[]
end

function status_GetTstartTstop(status::StepStatus)
    tstart, tstop = Ref{Cdouble}(), Ref{Cdouble}()
    @ccall libbraid.braid_StepStatusGetTstartTstop(status.ptr::Ptr{Cvoid}, tstart::Ref{Cdouble}, tstop::Ref{Cdouble})::Cint
    return tstart[], tstop[]
end

function status_GetTIndex(status::Union{StepStatus, AccessStatus})
    ti = Ref{Cint}()
    @ccall libbraid.braid_StatusGetTIndex(status.ptr::Ptr{Cvoid}, ti::Ref{Cint})::Cint
    return ti[]
end

function status_GetLevel(status::Union{StepStatus, AccessStatus})
    level = Ref{Cint}()
    @ccall libbraid.braid_StatusGetLevel(status.ptr::Ptr{Cvoid}, level::Ref{Cint})::Cint
    return level[]
end

function status_GetWrapperTest(status::Union{StepStatus, AccessStatus})
    wtest = Ref{Cint}()
    @ccall libbraid.braid_StatusGetWrapperTest(status.ptr::Ptr{Cvoid}, wtest::Ref{Cint})::Cint
    return (wtest[] > 0)
end

function status_GetIter(status::Union{StepStatus, AccessStatus})
    iter = Ref{Cint}()
    @ccall libbraid.braid_StatusGetIter(status.ptr::Ptr{Cvoid}, iter::Ref{Cint})::Cint
    return iter[]
end

function status_GetNLevels(status::Union{StepStatus, AccessStatus})
    nlevels = Ref{Cint}()
    @ccall libbraid.braid_StatusGetNLevels(status.ptr::Ptr{Cvoid}, nlevels::Ref{Cint})::Cint
    return nlevels[]
end

function status_GetNTPoints(status::Union{StepStatus, AccessStatus})
    ntpoints = Ref{Cint}()
    @ccall libbraid.braid_StatusGetNTPoints(status.ptr::Ptr{Cvoid}, ntpoints::Ref{Cint})::Cint
    return ntpoints[]
end

function status_GetResidual(status::Union{StepStatus, AccessStatus})
    res = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetResidual(status.ptr::Ptr{Cvoid}, res::Ref{Cdouble})::Cint
    return res[]
end

function status_GetDone(status::Union{StepStatus, AccessStatus})
    done = Ref{Cint}()
    @ccall libbraid.braid_StatusGetDone(status.ptr::Ptr{Cvoid}, done::Ref{Cint})::Cint
    return (done[] > 0)
end

function status_GetTILD(status::Union{StepStatus, AccessStatus})
    time = Ref{Cdouble}()
    iter = Ref{Cint}()
    level = Ref{Cint}()
    done = Ref{Cint}()
    @ccall libbraid.braid_StatusGetTILD(status.ptr::Ptr{Cvoid}, time::Ref{Cdouble}, iter::Ref{Cint}, level::Ref{Cint}, done::Ref{Cint})::Cint
    return time[], iter[], level[], (done[] > 0)
end

function status_GetCallingFunction(status::Union{StepStatus, AccessStatus})
    func = Ref{Cint}()
    @ccall libbraid.braid_StatusGetCallingFunction(status.ptr::Ptr{Cvoid}, func::Ref{Cint})::Cint
    return func[]
end

function status_GetTol(status::StepStatus)
    tol = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetTol(status.ptr::Ptr{Cvoid}, tol::Ref{Cdouble})::Cint
    return tol[]
end

function status_GetRNorms(status::StepStatus)
    iters = status_GetIter(status)
    nrequest = Ref{Cint}(iters)
    norms = zeros(Cdouble, iters)
    @ccall libbraid.braid_StatusGetRNorms(status.ptr::Ptr{Cvoid}, nrequest::Ref{Cint}, norms::Ref{Cdouble})::Cint
    return norms
end

function status_GetProc(status::Union{StepStatus, AccessStatus})
    proc = Ref{Cint}()
    @ccall libbraid.braid_StatusGetProc(status.ptr::Ptr{Cvoid}, proc::Ref{Cint})::Cint
    return proc[]
end

function status_GetDeltaRank(status::Union{StepStatus, AccessStatus})
    rank = Ref{Cint}()
    @ccall libbraid.braid_StatusGetDeltaRank(status.ptr::Ptr{Cvoid}, rank::Ref{Cint})::Cint
    return rank[]
end

# These are special cases
function status_GetLocalLyapExponents(status::AccessStatus)
    rank = status_GetDeltaRank(status)
    rank < 1 && return []

    exps = zeros(rank)
    num_retrieved = Ref(rank)
    @ccall  libbraid.braid_StatusGetLocalLyapExponents(status.ptr::Ptr{Cvoid}, exps::Ref{Cdouble}, num_retrieved::Ref{Cint})::Cint
    return exps
end

function status_GetBasisVectors(status::Union{StepStatus, AccessStatus})
    rank = status_GetDeltaRank(status)
    Ψ = []
    rank < 1 && return Ψ

    for i in 1:rank
        pp = malloc_null_double_ptr(Cvoid)
        GC.@preserve pp begin
            @ccall libbraid.braid_StatusGetBasisVec(status.ptr::Ptr{Cvoid}, pp::Ptr{Ptr{Cvoid}}, (i-1)::Cint)::Cint
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
