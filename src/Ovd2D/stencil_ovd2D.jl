get_functions(::Poisson) = (a_poisson, b_poisson, c_poisson, d_poisson, e_poisson)
get_functions(::FokkerPlanck) = (a_FP, b_FP, c_FP, d_FP, e_FP)

"""
    get_stencil(vars::Params, op_type::Operator, η::Float64)
    get_stencil(vars::Params, op_type::Operator)

Construct finite difference matrix for discretized operator (L or L^†)
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

    # get functions depending on operator choice
    a_fxn, b_fxn, c_fxn, d_fxn, e_fxn = get_functions(op_type)

	# main diagonal
    l = 1
    for i in 1:m, j in 1:m
        push!(is,l)
        push!(js,l)
        push!(vs, c_fxn(q[j], q[i], h, β, vars, η))
        l += 1
    end

	# upper and lower diagonals
    for i in 1:m, j in 1:m-1
        push!(is,(i-1)*m+j)
        push!(js,(i-1)*m+j+1)
        push!(vs, e_fxn(q[j], q[i], h, β, vars, η))

        push!(js,(i-1)*m+j)
        push!(is,(i-1)*m+j+1)
        push!(vs, a_fxn(q[j+1], q[i], h, β, vars, η))
    end

	# remaining 2 diagonals
    for i in 2:m, j in 1:m
        push!(is,(i-2)*m+j)
        push!(js,(i-1)*m+j)
        push!(vs, d_fxn(q[j], q[i-1], h, β, vars, η))

        push!(js,(i-2)*m+j)
        push!(is,(i-1)*m+j)
        push!(vs, b_fxn(q[j], q[i], h, β, vars, η))
    end

    # PBC terms
    for i in 1:m
        push!(is,(i-1)*m+1)
        push!(js,i*m)
        push!(vs, a_fxn(q[1], q[i], h, β, vars, η))

        push!(js,(i-1)*m+1)
        push!(is,i*m)
        push!(vs, e_fxn(q[end], q[i], h, β, vars, η))
    end

	# "corner" diagonal matrices in A
	for i in 1:m
        push!(is,i)
        push!(js,(m-1)*m+i)
        push!(vs, b_fxn(q[i], q[1], h, β, vars, η))

        push!(js,i)
        push!(is,(m-1)*m+i)
        push!(vs, d_fxn(q[i], q[end], h, β, vars, η))
	end

	return sparse(is, js, vs)
end

"Functions to build stencil for L (for Poisson)"
a_poisson(x, y, h, β::Float64, vars::Params, η::Float64=nothing) = -1/h^2/β - vars.d₁V(x,y)/2/h
b_poisson(x, y, h, β::Float64, vars::Params, η::Float64=nothing) = -1/h^2/β - vars.d₂V(x,y)/2/h
c_poisson(x, y, h, β::Float64, vars::Params, η::Float64=nothing) = 4/h^2/β
d_poisson(x, y, h, β::Float64, vars::Params, η::Float64=nothing) = -1/h^2/β + vars.d₂V(x,y)/2/h
e_poisson(x, y, h, β::Float64, vars::Params, η::Float64=nothing) = -1/h^2/β + vars.d₁V(x,y)/2/h

"Functions to build stencil for L^† (for Fokker-Planck)"
function a_FP(x, y, h, β, vars, η)
    # L part
    val = 1/h^2/β - vars.d₁V(x,y)/2/h
    # L̃ part
    val += η*vars.force(y)[1]/2/h

    return val
end

function b_FP(x, y, h, β, vars, η)
    # L part 
    val = 1/h^2/β - vars.d₂V(x,y)/2/h
    # L̃ part
    val += η*vars.force(y)[2]/2/h

    return val
end

function c_FP(x, y, h, β, vars, η)
    # L part 
    val = -4/h^2/β + vars.ΔV(x,y)
    # L̃ part (from nonconstant F)
    val += 0.0#-η*dforce(y)

    return val
end

function d_FP(x, y, h, β, vars, η)
    # L part 
    val = 1/h^2/β + vars.d₂V(x,y)/2/h
    # L̃ part
    val += -η*vars.force(y)[2]/2/h

    return val
end

function e_FP(x, y, h, β, vars, η)
    # L part 
    val = 1/h^2/β + vars.d₁V(x,y)/2/h
    # L̃ part
    val += -η*vars.force(y)[1]/2/h

    return val
end