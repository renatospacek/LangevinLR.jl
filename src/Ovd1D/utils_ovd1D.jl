# ===============================================================================
@with_kw mutable struct Params
    L::Float64
    m::Int64
    h::Float64
    q::Array{Float64,1}
    β::Float64
    V::Function
    dV::Function
    ΔV::Function
    force::Function
    gibbs::Function

    function Params(L, m, β, V, force)
        h = 2 * L / m
        q = range(-L, L - h + 1e-12, m)
        dV, ΔV = diff_potential(V)
        gibbs = _get_gibbs(V, L, β)

        new(L, m, h, q, β, V, dV, ΔV, force, gibbs)
    end
end

"Compute partial derivatives and Laplacian for a given potential `V`"
function diff_potential(V::Function)
    dV(x) = ForwardDiff.derivative(x -> V(x), x)
    ΔV(x) = ForwardDiff.derivative(x -> dV(x), x)

    return dV, ΔV
end

"Defines method for Gibbs distribution"
function _get_gibbs(V::Function, L::Real, β::Float64)
    Z, _ = quadgk(x -> exp(-β * V(x)), -L, L)

    return (x -> exp(-β * V(x)) / Z)
end

"Solve Ax = 0 for sparse A by finding eigenvec associated with the 0 eigenvalue"
function spectral_solver(A, vars::Params; normalize=true)
    # @unpack_Params vars

    val, vec = eigs(A, nev=1, which=:SM, tol=0.0)#, maxiter=2000)
    err = sqrt(vars.h) * norm(A * real.(vec) - val[1] * real.(vec))
    ψsm = real.(vec)
    normalize && ψsm ./= vars.h * sum(ψsm)

    return ψsm, err
end

"Solve linear system Ax = b with reshaped output"
function linear_solver(A, b)
    x = A \ b
    err = norm(A * x - b)
    # x, ch = lsqr(A, b, log=true); err = ch[:resnorm][end]

    return x, err
end