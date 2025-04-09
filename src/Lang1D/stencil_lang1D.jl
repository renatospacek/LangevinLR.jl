get_functions(::Poisson) = (a_poisson, b_poisson, c_poisson, d_poisson, e_poisson)
get_functions(::FokkerPlanck) = (a_FP, b_FP, c_FP, d_FP, e_FP)

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

    iV = Int[]
    jV = Int[]
    valV = Float64[]

    a_fxn, b_fxn, c_fxn, d_fxn, e_fxn = get_functions(op_type)

    # main diagonal
    l = 1
    for i in 1:mp, j in 1:mq
        push!(iV, l)
        push!(jV, l)
        push!(valV, c_fxn(q[j], p[i], hq, hp, γ, β, vars, η))
        l += 1
    end

    # upper and lower diagonals
    for i in 1:mp, j in 1:mq-1
        push!(iV, (i - 1) * mq + j)
        push!(jV, (i - 1) * mq + j + 1)
        push!(valV, e_fxn(q[j+1], p[i], hq, hp, γ, β, vars, η))

        push!(jV, (i - 1) * mq + j)
        push!(iV, (i - 1) * mq + j + 1)
        push!(valV, a_fxn(q[j], p[i], hq, hp, γ, β, vars, η))
    end

    # remaining 2 diagonals
    for i in 2:mp, j in 1:mq
        push!(iV, (i - 2) * mq + j)
        push!(jV, (i - 1) * mq + j)
        push!(valV, d_fxn(q[j], p[i-1], hq, hp, γ, β, vars, η))

        push!(jV, (i - 2) * mq + j)
        push!(iV, (i - 1) * mq + j)
        push!(valV, b_fxn(q[j], p[i], hq, hp, γ, β, vars, η))
    end

    # PBCs for q
    for i in 1:mp
        push!(iV, (i - 1) * mq + 1)
        push!(jV, i * mq)
        push!(valV, a_fxn(q[end], p[i], hq, hp, γ, β, vars, η))

        push!(jV, (i - 1) * mq + 1)
        push!(iV, i * mq)
        push!(valV, e_fxn(q[1], p[i], hq, hp, γ, β, vars, η))
    end

    return sparse(iV, jV, valV)
end

"Functions to build stencil for L (for Poisson)"
a_poisson(q, p, hq, hp, γ, β, vars::Params, η::Float64=nothing) = -p / 2 / hq
b_poisson(q, p, hq, hp, γ, β, vars::Params, η::Float64=nothing) = p / 2 / hp * γ + 1 / hp^2 * γ / β + vars.dV(q) / 2 / hp
c_poisson(q, p, hq, hp, γ, β, vars::Params, η::Float64=nothing) = -2 / hp^2 * γ / β
d_poisson(q, p, hq, hp, γ, β, vars::Params, η::Float64=nothing) = -p / 2 / hp * γ + 1 / hp^2 * γ / β - vars.dV(q) / 2 / hp
e_poisson(q, p, hq, hp, γ, β, vars::Params, η::Float64=nothing) = p / 2 / hq

"Functions to build stencil for L^† (for Fokker-Planck)"
function a_FP(q, p, hq, hp, γ, β, vars::Params, η)
    # L part
    val = p / 2 / hq

    return val
end

function b_FP(q, p, hq, hp, γ, β, vars::Params, η)
    # L part
    val = -p / 2 / hp * γ + 1 / hp^2 * γ / β - vars.dV(q) / 2 / hp
    # L̃ part
    val += vars.force(q) * η / 2 / hp

    return val
end

function c_FP(q, p, hq, hp, γ, β, vars::Params, η)
    # L part
    val = -2 / hp^2 * γ / β + γ

    return val
end

function d_FP(q, p, hq, hp, γ, β, vars::Params, η)
    # L part
    val = p / 2 / hp * γ + 1 / hp^2 * γ / β + vars.dV(q) / 2 / hp
    # L̃ part
    val += -vars.force(q) * η / 2 / hp

    return val
end

function e_FP(q, p, hq, hp, γ, β, vars::Params, η)
    # L part
    val = -p / 2 / hq

    return val
end