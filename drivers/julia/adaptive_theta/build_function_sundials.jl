using Symbolics, Dates

# code generation for sundials (C)
struct SunTarget <: Symbolics.BuildTargets end

function sunliterals(expr)
    expr isa Real && return :(RCONST($(float(expr))))
    return expr
end

function sunoperators(expr)
    expr isa Expr || return expr
    for e in expr.args
        if e isa Expr
            sunoperators(e)
        end
    end
    for i in eachindex(expr.args)
        if expr.args[i] isa Union{Float32, Float64, Rational}
            # RCONST is a macro from sundials which constructs a real constant
            expr.args[i] = :(RCONST($(float(expr.args[i]))))
        end
    end
    # Introduce another factor 1 to prevent contraction of terms like "5 * t" to "5t" (not valid C code)
    if expr.head==:call && expr.args[1]==:* && length(expr.args)==3 && isa(expr.args[2], Real) && isa(expr.args[3], Symbol)
        push!(expr.args, 1)
    # Power operator does not exist in C, replace by multiplication or "pow"
    elseif expr.head==:call && expr.args[1]==:^
        @assert length(expr.args)==3 "Don't know how to handle ^ operation with <> 2 arguments"
        x = expr.args[2]
        n = expr.args[3]
        empty!(expr.args)
        # Replace by multiplication/division if
        #   x is a symbol and  n is a small integer
        #   x is a more complex expression and n is ±1
        #   n is exactly 0
        if (isa(n,Integer) && ((isa(x, Symbol) && abs(n) <= 3) || abs(n) <= 1)) || n==0
            if n >= 0
                append!(expr.args, [:*, fill(x, n)...])
                # fill up with factor 1 so this expr can still be a multiplication
                while length(expr.args) < 3
                    push!(expr.args, 1)
                end
            else # inverse of the above
                if n==-1
                    term = x
                else
                    term = :( ($(x)) ^ ($(-n)))
                    coperators(term)
                end
                append!(expr.args, [:/, 1., term])
            end
        #... otherwise use "pow" function
        else
            append!(expr.args, [:SUNpowerI, x, n])
        end
    # replace bare real constants by RCONST
    elseif expr.head==:call && (expr.args[1]==:* || expr.args[1]==:+ || expr.args[1]==:/)
        for i in eachindex(expr.args)
            if expr.args[i] isa Real
                expr.args[i] = :(RCONST($(float(expr.args[i]))))
            end
        end
    end
    expr
end

function Symbolics._build_function(target::SunTarget, ex::AbstractArray, args...;
                                   conv     = Symbolics.toexpr,
                                   header   = false,
                                   fname    = :diffeqf,
                                   lhsname  = :du,
                                   rhsnames = [Symbol("RHS$i") for i in 1:length(args)])
    @info "Building function _$fname for target $target"
    fname = Symbol('_' * string(fname))
    ex = hcat([row for row in eachrow(ex)]...)
    varnumbercache = Symbolics.buildvarnumbercache(args...)
    equations = Vector{String}()
    for col ∈ axes(ex, 2), row ∈ axes(ex, 1)
        lhs = string(lhsname, "[", (col-1) * size(ex,1) + row-1, "]")
        rhs = Symbolics.numbered_expr(ex[row, col].val, varnumbercache, args...;
                                        lhsname  = lhsname,
                                        rhsnames = rhsnames,
                                        offset = -1) |> sunoperators |> sunliterals |> string
        push!(equations, string(lhs, " = ", rhs, ";"))
    end

    argstrs = join(vcat("sunrealtype* $(lhsname)",[typeof(args[i])<:AbstractArray ? "const sunrealtype* $(rhsnames[i])" : "const sunrealtype $(rhsnames[i])" for i in 1:length(args)]),", ")

    head = """
    #include <sundials/sundials_types.h>
    #include <sundials/sundials_math.h>

    """
    body = "void $fname($(argstrs...))\n{$([string("\n  ", eqn) for eqn ∈ equations]...)\n}\n"
    if header
        return head*body
    else
        return body
    end
end
