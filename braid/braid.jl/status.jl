# TODO: figure out status wrapper
function status_GetT(status)
    t = Ref{Cdouble}()
    @ccall libbraid.braid_StatusGetT(status::Ptr{Cvoid}, t::Ref{Cdouble})::Cint
    return t[]
end

function status_GetTIndex(status)
    ti = Ref{Cint}()
    @ccall libbraid.braid_StatusGetTIndex(status::Ptr{Cvoid}, ti::Ref{Cint})::Cint
    return ti[]
end

# These are special cases
function status_GetLocalLyapExponents(status)
    rank = Ref{Cint}()
    @ccall libbraid.braid_StatusGetDeltaRank(status::Ptr{Cvoid}, rank::Ref{Cint})::Cint
    rank[] < 1 && return []

    exps = zeros(rank[])
    num_retrieved = Ref(rank[])
    @ccall  libbraid.braid_StatusGetLocalLyapExponents(status::Ptr{Cvoid}, exps::Ref{Cdouble}, num_retrieved::Ref{Cint})::Cint
    return exps
end

function status_GetBasisVectors(status)
    rank = Ref{Cint}()
    @ccall libbraid.braid_StatusGetDeltaRank(status::Ptr{Cvoid}, rank::Ref{Cint})::Cint
    Ψ = []
    rank[] < 1 && return Ψ

    for i in 1:rank[]
        # double pointer to NULL
        pp = get_null_double_ptr(Cvoid)
        @ccall libbraid.braid_StatusGetBasisVec(status::Ptr{Cvoid}, pp::Ptr{Ptr{Cvoid}}, (i-1)::Cint)::Cint
        p = unsafe_load(pp)
        if p !== C_NULL
            φ = unsafe_pointer_to_objref(p)
            push!(Ψ, φ.user_vector)
        end
    end

    return Ψ
end
