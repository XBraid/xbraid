using ForwardDiff, LinearAlgebra, PreallocationTools
using IterativeSolvers, LinearMaps, SparseArrays

"""
struct to package preallocated caches for storing temp coords and velocity
as well as the sparse poisson matrix needed by project_incompressible!()
"""
struct TGApp
    x_d::DiffCache
    u_d::DiffCache
    ϕ_d::DiffCache
    Δ_s::SparseMatrixCSC
    solution
    lyapunov_vecs
    lyapunov_exps
end
TGApp(x::AbstractArray, u::AbstractArray, ϕ::AbstractArray, Δ::SparseMatrixCSC) = TGApp(DiffCache(x), DiffCache(u), DiffCache(ϕ), Δ, [], [], [])

"""
reshapes 1D array into vector of 3D arrays, with no allocation
"""
function get_views(u)
    @views begin
        u_x = reshape(u[1:nₓ^3], (nₓ, nₓ, nₓ))
        u_y = reshape(u[nₓ^3+1:2*nₓ^3], (nₓ, nₓ, nₓ))
        u_z = reshape(u[2*nₓ^3+1:end], (nₓ, nₓ, nₓ))
    end
    return [u_x, u_y, u_z]
end

# some helpers for doing modulo arithmetic with CartesianIndex
Base.mod1(I::CartesianIndex, J::CartesianIndex) = CartesianIndex(broadcast(mod1, Tuple(I), Tuple(J)))
δ(I::CartesianIndex{N}, dim) where N = CartesianIndex(ntuple(i -> i == dim ? 1 : 0, N))

"""
differentiable trilinear interpolation used by semi-lagrangian discretization
"""
function interp_linear_periodic(x, y, z, F)
    bound = CartesianIndex(nₓ, nₓ, nₓ)
    # scale coordinates
    xs, ys, zs = x/Δx, y/Δx, z/Δx
    # get integer parts and remainders
    # that is, the origin vertex of the current grid cell,
    # and the relative coordinates from that vertex
    i, j, k = trunc.((xs, ys, zs))
    r, s, t = xs - i, ys - j, zs - k
    i, j, k = Int.([i, j, k] .+ 1)

    I = mod1(CartesianIndex(i, j, k), bound)
    I_ones = CartesianIndex(1, 1, 1)
    C = zeros(eltype(F), 2, 2, 2)
    for I_c in CartesianIndices(C)
        I_rel = I_c - I_ones
        C[I_c] = F[mod1(I + I_rel, bound)]
    end

    # interp over x, y, then z
    C = (1-r)*C[1,:,:] + r*C[2,:,:]
    C = (1-s)*C[1,:]   + s*C[2,:]
    C = (1-t)*C[1]     + t*C[2]
    return C
end

function TaylorGreen!(u)
    u_x, u_y, u_z = u
    x, y, z = coords
    @. u_x =  sin(x) * cos(y) * cos(z)
    @. u_y = -cos(x) * sin(y) * cos(z)
    @. u_z = 0.
    return
end

"""
finite difference approx first derivative
of A wrt dimension *dim*
"""
function ∂(A, dim;h=Δx)
    # This works for A of arbitrary dimension
    out = similar(A)
    R = CartesianIndices(A)
    bound = last(R)
    for I ∈ R
        I⁻ = mod1(I - δ(I, dim), bound)
        # second order
        I⁺ = mod1(I + δ(I, dim), bound)
        out[I] = 1/2h * (A[I⁺] - A[I⁻])
        # first order
        # out[I] = 1/h * (A[I] - A[I⁻])
    end
    return out
end

# divergence, curl, gradient, and laplacian
∇_dot(V) = sum(enumerate(V)) do (i, V_i)
    ∂(V_i, i)
end
function ∇X(V)
    out_arr = zeros(3*nₓ^3)
    out = get_views(out_arr)
    Vx, Vy, Vz = V
    x, y, z = 1:3
    out[x] .= ∂(Vz, y) - ∂(Vy, z)
    out[y] .= ∂(Vx, z) - ∂(Vz, x)
    out[z] .= ∂(Vy, x) - ∂(Vx, y)
    return out
end
∇(F) = [∂(F, i) for i ∈ 1:3]
# Δ(F) = ∇_dot(∇(F)) # this is asymmetric...


function my_poisson(f_arr)
    out_arr = similar(f_arr)
    F = reshape(f_arr, (nₓ, nₓ, nₓ))
    out = reshape(out_arr, (nₓ, nₓ, nₓ))
    R = CartesianIndices(F)
    bound = last(R)
    @inbounds for I ∈ R
        out[I] = -6*F[I]
        @inbounds for dim ∈ 1:3
            I⁺ = mod1(I + δ(I, dim), bound)
            I⁻ = mod1(I - δ(I, dim), bound)
            out[I] += F[I⁺] + F[I⁻]
        end
        out[I] *= 1/Δx^2
    end
    # need a single point condition u(0,0,0) = 0., else
    # this system is singular (ones(nₓ^3) nullspace)
    out[first(R)] += F[first(R)]
    return out_arr
end

function my_poisson_mat(nₓ)
    flat_index(I::CartesianIndex) = I[1] + (I[2]-1)*nₓ + (I[3]-1)*nₓ^2
    rowinds = Array{Int}([])
    colinds = Array{Int}([])
    nzs = Array{Float}([])
    function make_connection!(I_r, I_c, nz)
        push!(rowinds, flat_index(I_r))
        push!(colinds, flat_index(I_c))
        push!(nzs, nz)
    end
    diag = -6/Δx^2
    offd = 1/Δx^2
    bound = CartesianIndex((nₓ, nₓ, nₓ))
    for I in CartesianIndices((1:nₓ, 1:nₓ, 1:nₓ))
        make_connection!(I, I, diag)
        for dim ∈ 1:3
            I⁺ = mod1(I + δ(I, dim), bound)
            I⁻ = mod1(I - δ(I, dim), bound)
            make_connection!(I, I⁺, offd)
            make_connection!(I, I⁻, offd)
        end
    end
    # point condition u(0,0,0) = 0., else system is singular
    nzs[1] += 1.
    return sparse(rowinds, colinds, nzs)
end

"""
compute self advection of u
"""
function advect_semi_lagrangian!(app::TGApp, u_arr, Δt)
    # get cached arrays (either real or dual, depending on typeof(u))
    coords_new_arr = get_tmp(app.x_d, u_arr)
    u_new_arr = get_tmp(app.u_d, u_arr)

    u = get_views(u_arr)
    u_new = get_views(u_new_arr)

    # (1) trace each grid point backwards along u
    @. coords_new_arr = coords_arr - Δt * u_arr
    coords_new = get_views(coords_new_arr)

    # (2) interpolate the velocity field at the new grid points
    for i in 1:3
        u_new[i] .= interp_linear_periodic.(coords_new..., Ref(u[i]))
    end
    u_arr .= u_new_arr
    return u_arr
end

"""
Solve Poisson problem for incompressible velocity field
Δϕ = ∇⋅u
u = u - ∇ϕ
"""
function project_incompressible!(app::TGApp, u_arr)
    u = get_views(u_arr)
    rhs = reshape(∇_dot(u), nₓ^3)
    ϕ_arr = cg(Δ_lin, rhs)
    ϕ = reshape(ϕ_arr, (nₓ, nₓ, nₓ))
    for i in 1:3
        u[i] .-= ∂(ϕ, i)
    end
    return u_arr
end

function step(app::TGApp, u_arr, Δt)
    advect_semi_lagrangian!(u_arr, Δt)
    # diffuse_backward_Euler!(u_arr, Δt) # numerical diffusion may be enough
    project_incompressible!(u_arr)
    return u_arr
end