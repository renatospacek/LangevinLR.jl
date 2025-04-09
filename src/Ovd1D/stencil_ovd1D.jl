get_functions(::Poisson) = (a_poisson, b_poisson, c_poisson)
get_functions(::FokkerPlanck) = (a_FP, b_FP, c_FP)

"""
    get_stencil(vars::Params, op_type::Operator, η::Float64)
    get_stencil(vars::Params, op_type::Operator)

Constructs finite difference matrix for discretized operator (L or L^†)
    - for `op_type::Poisson`, returns generator L
    - for `op_type::FokkerPlanck`, returns L² adjoint L^†

η is set to 0 if not passed as an argument
"""
get_stencil(vars::Params, op_type::Operator) = get_stencil(vars::Params, op_type::Operator, 0.0)

function get_stencil(vars::Params, op_type::Operator, η::Float64)
    @unpack_Params vars

    is = Int[]
    js = Int[]
    vs = Float64[]

    a_fxn, b_fxn, c_fxn = get_functions(op_type)

    # main diagonal
    for i in 1:m
        push!(is, i)
        push!(js, i)
        push!(vs, b_fxn(q[i], h, β, vars, η))
    end

    # off diagonals
    for i in 1:m-1
        push!(is, i)
        push!(js, i + 1)
        push!(vs, c_fxn(q[i], h, β, vars, η))

        push!(is, i + 1)
        push!(js, i)
        push!(vs, a_fxn(q[i+1], h, β, vars, η))
    end

    # PBCs
    push!(is, 1)
    push!(js, m)
    push!(vs, a_fxn(q[1], h, β, vars, η))

    push!(is, m)
    push!(js, 1)
    push!(vs, c_fxn(q[m], h, β, vars, η))

    return sparse(is, js, vs)
end

"Functions to build stencil for L (for Poisson)"
a_poisson(q, h, β, vars::Params, η::Float64=nothing) = -1 / h^2 / β - vars.dV(q) / 2 / h
b_poisson(q, h, β, vars::Params, η::Float64=nothing) = 2 / h^2 / β
c_poisson(q, h, β, vars::Params, η::Float64=nothing) = -1 / h^2 / β + vars.dV(q) / 2 / h

"Functions to build stencil for L^† (for Fokker-Planck)"
function a_FP(q, h, β, vars::Params, η)
    # L_0 part
    val = 1 / h^2 - vars.dV(q) / 2 / h
    # L̃ part
    val += vars.force(q) * η / 2 / h

    return val
end

function b_FP(q, h, β, vars::Params, η)
    # L_0 part
    val = -2 / h^2 + vars.ΔV(q)

    return val
end

function c_FP(q, h, β, vars::Params, η)
    # L_0 part
    val = 1 / h^2 + vars.dV(q) / 2 / h
    # L̃ part
    val += -vars.force(q) * η / 2 / h

    return val
end