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


module Status
export AccessStatus, StepStatus, SyncStatus
export CallerError, CallerFInterp, CallerFRestrict, CallerFRefine, CallerFAccess
export CallerFRefine_AfterInitHier, CallerDrive_TopCycle, CallerFCRelax
export CallerDrive_AfterInit, CallerBaseStep_diff, CallerComputeFullRNorm
export CallerFASResidual, CallerResidual, CallerInitGuess, CallerWarmup

using ..XBraid.BraidUtils: libbraid, malloc_null_double_ptr
using ..XBraid: BraidVector

# using MPI: MPI_Comm, Comm

# calling functions
const CallerError = -1
const CallerFInterp = 0
const CallerFRestrict = 1
const CallerFRefine = 2
const CallerFAccess = 3
const CallerFRefine_AfterInitHier = 4
const CallerDrive_TopCycle = 5
const CallerFCRelax = 6
const CallerDrive_AfterInit = 7
const CallerBaseStep_diff = 8
const CallerComputeFullRNorm = 9
const CallerFASResidual = 10
const CallerResidual = 11
const CallerInitGuess = 12
const CallerWarmup = 13

# status structures
struct StepStatus
    ptr::Ptr{Cvoid}
end
StepStatus() = StepStatus(C_NULL)

struct AccessStatus
    ptr::Ptr{Cvoid}
end
AccessStatus() = AccessStatus(C_NULL)

struct SyncStatus
    ptr::Ptr{Cvoid}
end
SyncStatus() = SyncStatus(C_NULL)

# status get/set functions (not exported)
# @ccall requires these type annotations to be known at compile time, so macros won't work here
function getT(status::Union{StepStatus, AccessStatus})
    status.ptr == C_NULL && return 0.0
    t = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetT(status.ptr::Ptr{Cvoid}, t::Ref{Cdouble})::Cint
    return t[]
end

function getTStop(status::StepStatus)
    status.ptr == C_NULL && return 0.0
    tstop = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetTStop(status.ptr::Ptr{Cvoid}, tstop::Ref{Cdouble})::Cint
    return tstop[]
end

function getTstartTstop(status::StepStatus)
    status.ptr == C_NULL && return (0.0, 0.0)
    tstart, tstop = Ref{Cdouble}(), Ref{Cdouble}()
    @ccall libbraid.braid_StepStatusGetTstartTstop(status.ptr::Ptr{Cvoid}, tstart::Ref{Cdouble}, tstop::Ref{Cdouble})::Cint
    return tstart[], tstop[]
end

function getTIndex(status::Union{StepStatus, AccessStatus})
    status.ptr == C_NULL && return 0
    ti = Ref{Cint}()
    @ccall libbraid.braid_StatusGetTIndex(status.ptr::Ptr{Cvoid}, ti::Ref{Cint})::Cint
    # julia uses 1-based indexing
    return ti[] + 1
end

function getLevel(status::Union{StepStatus, AccessStatus, SyncStatus})
    status.ptr == C_NULL && return 0
    level = Ref{Cint}()
    @ccall libbraid.braid_StatusGetLevel(status.ptr::Ptr{Cvoid}, level::Ref{Cint})::Cint
    return level[]
end

function getCFactor(status::Union{StepStatus, AccessStatus}, level::Integer)
    status.ptr == C_NULL && return 1
    cfactor = Ref{Cint}()
    @ccall libbraid.braid_StatusGetCFactor(status.ptr::Ptr{Cvoid}, cfactor::Ref{Cint}, level::Cint)::Cint
    return cfactor[]
end

function getCFactor(status::Union{StepStatus, AccessStatus})
    status.ptr == C_NULL && return 0
    getCFactor(status, getLevel(status))
end

function getWrapperTest(status::Union{StepStatus, AccessStatus})
    status.ptr == C_NULL && return false
    wtest = Ref{Cint}()
    @ccall libbraid.braid_StatusGetWrapperTest(status.ptr::Ptr{Cvoid}, wtest::Ref{Cint})::Cint
    return (wtest[] > 0)
end

function getIter(status::Union{StepStatus, AccessStatus, SyncStatus})
    status.ptr == C_NULL && return 0
    iter = Ref{Cint}()
    @ccall libbraid.braid_StatusGetIter(status.ptr::Ptr{Cvoid}, iter::Ref{Cint})::Cint
    return iter[]
end

function getNLevels(status::Union{StepStatus, AccessStatus, SyncStatus})
    status.ptr == C_NULL && return 0
    nlevels = Ref{Cint}()
    @ccall libbraid.braid_StatusGetNLevels(status.ptr::Ptr{Cvoid}, nlevels::Ref{Cint})::Cint
    return nlevels[]
end

function setRFactor(status::StepStatus, rfactor::Real)
    status.ptr == C_NULL && return nothing
    @ccall libbraid.braid_StatusSetRFactor(status.ptr::Ptr{Cvoid}, rfactor::Cdouble)::Cint
    return nothing
end

function getNRefine(status::Union{StepStatus, AccessStatus, SyncStatus})
    status.ptr == C_NULL && return 0
    nrefine = Ref{Cint}()
    @ccall libbraid.braid_StatusGetNRefine(status.ptr::Ptr{Cvoid}, nrefine::Ref{Cint})::Cint
    return nrefine[]
end

function getNTPoints(status::Union{StepStatus, AccessStatus, SyncStatus})
    status.ptr == C_NULL && return 0
    ntpoints = Ref{Cint}()
    @ccall libbraid.braid_StatusGetNTPoints(status.ptr::Ptr{Cvoid}, ntpoints::Ref{Cint})::Cint
    return ntpoints[]
end

function getResidual(status::Union{StepStatus, AccessStatus})
    status.ptr == C_NULL && return 0.0
    res = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetResidual(status.ptr::Ptr{Cvoid}, res::Ref{Cdouble})::Cint
    return res[]
end

function getDone(status::Union{StepStatus, AccessStatus, SyncStatus})
    status.ptr == C_NULL && return false
    done = Ref{Cint}()
    @ccall libbraid.braid_StatusGetDone(status.ptr::Ptr{Cvoid}, done::Ref{Cint})::Cint
    return (done[] > 0)
end

function getTILD(status::Union{StepStatus, AccessStatus})
    status.ptr == C_NULL && return (0.0, 0, 0, false)
    time = Ref{Cdouble}()
    iter = Ref{Cint}()
    level = Ref{Cint}()
    done = Ref{Cint}()
    @ccall libbraid.braid_StatusGetTILD(status.ptr::Ptr{Cvoid}, time::Ref{Cdouble}, iter::Ref{Cint}, level::Ref{Cint}, done::Ref{Cint})::Cint
    return time[], iter[], level[], (done[] > 0)
end

function getTIUL(status::Union{StepStatus, SyncStatus}, level::Integer)
    status.ptr == C_NULL && return (0, 0)
    iu = Ref{Cint}()
    il = Ref{Cint}()
    @ccall libbraid.braid_StatusGetTIUL(status.ptr::Ptr{Cvoid}, iu::Ref{Cint}, il::Ref{Cint}, level::Cint)::Cint
    # julia uses 1-based indexing
    return iu[] + 1, il[] + 1
end

function getCallingFunction(status::Union{StepStatus, AccessStatus, SyncStatus})
    status.ptr == C_NULL && return CallerError
    func = Ref{Cint}()
    @ccall libbraid.braid_StatusGetCallingFunction(status.ptr::Ptr{Cvoid}, func::Ref{Cint})::Cint
    return func[]
end

function getTol(status::StepStatus)
    status.ptr == C_NULL && return 0.0
    tol = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetTol(status.ptr::Ptr{Cvoid}, tol::Ref{Cdouble})::Cint
    return tol[]
end

function getRNorms(status::StepStatus)
    status.ptr == C_NULL && return zeros(0)
    iters = getIter(status)
    nrequest = Ref{Cint}(iters)
    norms = zeros(Cdouble, iters)
    @ccall libbraid.braid_StatusGetRNorms(status.ptr::Ptr{Cvoid}, nrequest::Ref{Cint}, norms::Ref{Cdouble})::Cint
    return norms
end

function getProc(status::Union{StepStatus, AccessStatus, SyncStatus})
    status.ptr == C_NULL && return 0
    proc = Ref{Cint}()
    @ccall libbraid.braid_StatusGetProc(status.ptr::Ptr{Cvoid}, proc::Ref{Cint})::Cint
    return proc[]
end

# function getTComm(status::SyncStatus)
#     status.ptr == C_NULL && return MPI.Comm(MPI.COMM_NULL)
#     tcomm = Ref{MPI.MPI_Comm}()
#     @ccall libbraid.braid_StatusGetTComm(status.ptr::Ptr{Cvoid}, tcomm::Ref{MPI.MPI_Comm})::Cint
#     return MPI.Comm(tcomm[])
# end

function getDeltaRank(status::Union{StepStatus, AccessStatus})
    status.ptr == C_NULL && return 0
    rank = Ref{Cint}()
    @ccall libbraid.braid_StatusGetDeltaRank(status.ptr::Ptr{Cvoid}, rank::Ref{Cint})::Cint
    return rank[]
end

# These are special cases
function getLocalLyapExponents(status::AccessStatus)
    status.ptr == C_NULL && return zeros(0)
    rank = getDeltaRank(status)
    rank < 1 && return []

    exps = zeros(rank)
    num_retrieved = Ref(rank)
    @ccall  libbraid.braid_StatusGetLocalLyapExponents(status.ptr::Ptr{Cvoid}, exps::Ref{Cdouble}, num_retrieved::Ref{Cint})::Cint
    return exps
end

function getBasisVectors(status::Union{StepStatus, AccessStatus})
    status.ptr == C_NULL && return []
    rank = getDeltaRank(status)
    rank < 1 && return []

    Ψ = []
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

end # module Status