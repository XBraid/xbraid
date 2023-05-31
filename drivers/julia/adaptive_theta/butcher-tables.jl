using DataStructures, PrettyTables
using Symbolics, NLsolve

include("build_function_sundials.jl")

struct ButcherTable
    A::AbstractMatrix
    b::AbstractVector
    c::AbstractVector
    s::Int
end
ButcherTable(A, b) = ButcherTable(A, b, A*ones(length(b)), length(b))
ButcherTable(A, b, c) = ButcherTable(A, b, c, length(b))
function ButcherTable(A::Vector{<:Real}, b::Vector{<:Real}, c::Vector{<:Real}, s::Integer)
    return ButcherTable(reshape(A, s, s), b, c, s)
end
function Base.show(io::IO, bt::ButcherTable)
    print(io, "ButcherTable: \n")
    # pretty table
    b = [1.0; bt.b]
    data = [bt.c bt.A; b']
    pretty_table(
        data;
        show_header=false,
        body_hlines=[bt.s],
        vlines=[:begin, 1, :end]
    )
end

forest = OrderedDict(
    # 2nd order
    :{τ} => 1,
    # 3rd order
    :{{τ}} => 2,
    :{τ, τ} => 3,
    # 4th order
    :{{{τ}}} => 4,
    :{{τ, τ}} => 5,
    :{τ, {τ}} => 6,
    :{τ, τ, τ} => 7,
    # 5th order
    :{{{{τ}}}} => 8,
    :{{{τ, τ}}} => 9,
    :{{τ, {τ}}} => 10,
    :{{τ, τ, τ}} => 11,
    :{τ, {{τ}}} => 12,
    :{τ, {τ, τ}} => 13,
    :{{τ}, {τ}} => 14,
    :{τ, τ, {τ}} => 15,
    :{τ, τ, τ, τ} => 16
)

elem_weights = [
    # 2nd order
    (A, b, c) -> b' * c,
    # 3rd order
    (A, b, c) -> b' * A * c,
    (A, b, c) -> b' * c.^2,
    # 4th order
    (A, b, c) -> b' * A * A * c,
    (A, b, c) -> b' * A * c.^2,
    (A, b, c) -> b' * (c .* (A*c)),
    (A, b, c) -> b' * c.^3,
    # 5th order
    (A, b, c) -> b' * A * A * A * c,
    (A, b, c) -> b' * A * A * c.^2,
    (A, b, c) -> b' * A * (c .* (A*c)),
    (A, b, c) -> b' * A * c.^3,
    (A, b, c) -> b' * (c .* (A*A*c)),
    (A, b, c) -> b' * (c .* (A*c.^2)),
    (A, b, c) -> b' * (A * c).^2,
    (A, b, c) -> b' * (c .* (c .* (A*c))),
    (A, b, c) -> b' * c.^4
]

rhs_classical = [
    # 2nd order
    1/2,
    # 3rd order
    1/6, 1/3,
    # 4th order
    1/24, 1/12, 1/8, 1/4,
    # 5th order
    1/120, 1/60, 1/40, 1/20, 1/30, 1/15, 1/20, 1/10, 1/5
]
num_conditions = [0, 1, 3, 7, 16]

ψ(B::ButcherTable, tree::Expr) = elem_weights[forest[tree]](B.A, B.b, B.c)

@variables θ[1:15]

function gen_lhs_func(B::ButcherTable, order::Integer)
    num_conds = num_conditions[order]
    fA, fb, fc = [eval(build_function(syms, collect(θ[1:num_conds]); checkbounds=true)[1]) for syms ∈ (B.A, B.b, B.c)]
    conds_sym = [expand(ψ(B, t)) for t ∈ forest.keys[1:num_conds]]
    conds_jac = Symbolics.jacobian(conds_sym, θ[1:num_conds])
    f = eval(build_function(conds_sym, collect(θ[1:num_conds]); checkbounds=true)[2])
    J = eval(build_function(conds_jac, collect(θ[1:num_conds]); checkbounds=true)[2])
    function fill(θ::AbstractVector)
        @assert length(θ) == num_conds
        ButcherTable(fA(θ), fb(θ), fc(θ), B.s)
    end
    return f, J, fill
end

function solve_order_conditions(f!::Function, J!::Function, fill::Function, order::Integer, rhs::Vector{<:Real};
                                guess=zeros(num_conditions[order]), method=:trust_region)
    @assert length(rhs) == num_conditions[order]
    function g!(G, θ)
        f!(G, θ)
        G .-= rhs
    end
    result = nlsolve(g!, J!, guess; method=method)
    if !result.f_converged
        @warn "Order conditions solver not converged"
        @info result
        result = nlsolve(g!, J!, guess; method=:newton)
    end
    fill(result.zero)
end

function gen_c_lhs_func(B::ButcherTable, order::Integer, name::AbstractString)
    rhsnames = [:th]
    num_conds = num_conditions[order]
    # fill butcher table
    fA = build_function(
        B.A, collect(θ[1:num_conds]);
        target=SunTarget(),
        fname=name * "_btable_A",
        lhsname=:A, rhsnames=rhsnames
    )
    fb = build_function(
        B.b, collect(θ[1:num_conds]);
        target=SunTarget(),
        fname=name * "_btable_b",
        lhsname=:b, rhsnames=rhsnames
    )
    fc = build_function(
        B.c, collect(θ[1:num_conds]);
        target=SunTarget(),
        fname=name * "_btable_c",
        lhsname=:c, rhsnames=rhsnames
    )
    # compute order conditions
    conds_sym = [expand(ψ(B, t)) for t ∈ forest.keys[1:num_conds]]
    conds_jac = Symbolics.jacobian(conds_sym, θ[1:num_conds])
    func = build_function(
        conds_sym, collect(θ[1:num_conds]);
        target=SunTarget(), header=true,
        fname=name * "_lhs",
        lhsname=:phi, rhsnames=rhsnames
    )
    func_j = build_function(
        conds_jac, collect(θ[1:num_conds]);
        target=SunTarget(),
        fname=name * "_lhs_jac",
        lhsname=:phi_J, rhsnames=rhsnames
    )
    open("c_funcs/$name.c", "w") do file
        println(file, func)
        println(file, func_j)
        println(file, fA)
        println(file, fb)
        println(file, fc)
    end
end

beuler = ButcherTable([1.0], [1.0], [1.0])
sdirk212 = ButcherTable([1. 0.; -1. 1.], [1 / 2, 1 / 2], [1., 0.])
sdirk212_emb = [1., 0.]

#sdirk33
x = 0.4358665215
sdirk33 = ButcherTable([x 0.0 0.0; (1-x)/2 x 0; (-3.0x^2/2 + 4.0x - 1/4) (3x^2/2 - 5.0x + 5/4) x], [(-3.0x^2/2 + 4.0x - 1/4), (3.0x^2/2 - 5.0x + 5/4), x])

#esdirk423
diag = 1767732205903 / 4055673282236
A3 = zeros(4, 4)
A3[2, 1] = diag
A3[2, 2] = diag
A3[3, 1] = 2746238789719 / 10658868560708
A3[3, 2] = -640167445237 / 6845629431997
A3[3, 3] = diag
A3[4, 1] = 1471266399579 / 7840856788654
A3[4, 2] = -4482444167858 / 7529755066697
A3[4, 3] = 11266239266428 / 11593286722821
A3[4, 4] = diag
esdirk3 = ButcherTable(A3, A3[4, :], [0, 2diag, 3 / 5, 1])
esdirk3_emb = [2756255671327/12835298489170, -10771552573575/22201958757719, 9247589265047/10645013368117, 2193209047091/5459859503100]

#sdirk534
A4 = [
        1/4       0        0       0     0;
        1/2       1/4      0       0     0;
        17/50    -1/25     1/4     0     0;
        371/1360 -137/2720 15/544  1/4   0;
        25/24    -49/48    125/16 -85/12 1/4
]
sdirk4 = ButcherTable(A4, A4[5, :])


# θ methods

# can approximate up to 2nd order
θesdirk2 = ButcherTable([0.0 0.0; 1-θ[1] θ[1]], [1-θ[1], θ[1]])
θesdirk2_lhs!, θesdirk2_J!, θesdirk2_fill = gen_lhs_func(θesdirk2, 2)
θesdirk2_guess = [1/2]
gen_c_lhs_func(θesdirk2, 2, "theta_esdirk2")
#=
    0 │     0   0 
    1 │ 1 - θ   θ 
──────┼───────────
    1 │ 1 - θ   θ 
=#

# can approximate up to 2nd order
θsdirk2 = ButcherTable([θ[1] 0.0; 1-θ[1] θ[1]], [1-θ[1], θ[1]])
θsdirk2_lhs!, θsdirk2_J!, θsdirk2_fill = gen_lhs_func(θsdirk2, 2)
θsdirk2_guess = [1-√2/2]
gen_c_lhs_func(θsdirk2, 2, "theta_sdirk2")
#=
    θ │     θ   0 
    1 │ 1 - θ   θ 
──────┼───────────
    1 │ 1 - θ   θ 
=#

# can approximate up to 3rd order
θesdirk3 = ButcherTable([0. 0. 0.; θ[2]-θ[1] θ[1] 0; θ[3] (1.0-θ[3]-θ[1]) θ[1]], [θ[3], 1.0-θ[3]-θ[1], θ[1]])
θesdirk3_lhs!, θesdirk3_J!, θesdirk3_fill = gen_lhs_func(θesdirk3, 3)
gen_c_lhs_func(θesdirk3, 3, "theta_esdirk3")
#=
   0  │      0            0   0
   θ₂ │ θ₂ - θ₁           θ₁  0
   1  │      θ₃  1 - θ₃ - θ₁  θ₁
──────┼──────────────────────────
   2  │      θ₃  1 - θ₃ - θ₁  θ₁
=#

θsdirk3 = ButcherTable([θ[1] 0.0 0.0; θ[2]-θ[1] θ[1] 0; (1.0-θ[3]-θ[1]) θ[3] θ[1]], [1.0-θ[3]-θ[1], θ[3], θ[1]])
θsdirk3_lhs!, θsdirk3_J!, θsdirk3_fill = gen_lhs_func(θsdirk3, 3)
θsdirk3_guess = [0.4358665215, 0.71793326075, -0.6443631706532353]
gen_c_lhs_func(θsdirk3, 3, "theta_sdirk3")
#=
   θ₁ │          θ₁  0   0
   θ₂ │     θ₂ - θ₁  θ₁  0
   1  │ 1 - θ₃ - θ₁  θ₃  θ₁
──────┼──────────────────────────
   2  │ 1 - θ₃ - θ₁  θ₃  θ₁
=#

# can approximate up to 4th order
# θsdirk4 = ButcherTable([θ[1] 0 0 0 0; θ[2]-θ[1] θ[1] 0 0; θ[3] (1-θ[3]-θ[1]) θ[1] 0 0; θ[4] 0 (1-θ[4]-θ[3]-θ[1]) θ[1]], [θ[4], 0, 1-θ[4]-θ[3]-θ[1], θ[1]])
# θsdirk4_lhs, θsdirk4_J = gen_lhs_func(θsdirk4, 4)
# gen_c_lhs_func(θsdirk4, 4, "theta_sdirk4")
# #=
#    θ₁ │      θ₁           0   0   0
#    θ₂ │ θ₂ - θ₁           θ₁  0   0
#    θ₃ │      θ₃  1 - θ₃ - θ₁  θ₁  0
#    1  │      θ₄           0   θ₃  θ₁
