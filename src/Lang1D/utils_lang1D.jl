# ===============================================================================
@with_kw mutable struct Params
    Lq::AbstractFloat
    Lp::Float64
    mq::Int64
    mp::Int64
    hq::Float64
    hp::Float64
    q::Array{Float64,1}
    p::Array{Float64,1}
    X::Array{Float64,2}
    Y::Array{Float64,2}
    γ::Float64
    β::Float64
    V::Function
    dV::Function
    ΔV::Function
    force::Function
    gibbs::Function

    function Params(Lp, mq::Int, mp::Int, γ, β, V, force)
        Lq = π
        hq = 2Lq / mq
        hp = 2Lp / (mp - 1)

        q = range(-Lq, Lq - hq + 1e-12, mq)
        p = range(-Lp, Lp, mp)
        X0, Y0 = meshgrid(q, p)

        dV, ΔV = diff_potential(V)
        gibbs = _get_gibbs(V, Lq, β)

        X = X0'
        Y = Y0'

        new(Lq, Lp, mq, mp, hq, hp, q, p, X, Y, γ, β, V, dV, ΔV, force, gibbs)
    end
end

"Compute partial derivatives and Laplacian for a given potential `V`"
function diff_potential(V::Function)
    dV(x) = ForwardDiff.derivative(x -> V(x), x)
    ΔV(x) = ForwardDiff.derivative(x -> dV(x), x)

    return dV, ΔV
end

"Defines method for Gibbs distribution"
function _get_gibbs(V::Function, Lq::Real, β::Float64)
    # integral_p = sqrt(2π / β)
    # integral_q, _ = quadgk(x -> exp(-β * V(x)), -π, π)
    # Z = integral_p * integral_q
    Z = hcubature(x -> exp(-β * V(x[1]) - β * x[2]^2 / 2), [-Lq, -6.0], [Lq, 6.0])[1]

    return ((x, y) -> exp(-β * V(x) - β * y^2 / 2) / Z)
end

"Solve Ax = 0 for sparse A by finding eigenvec associated with the 0 eigenvalue"
function spectral_solver(A, vars::Params; normalize=true)
    # @unpack_Params vars

    val, vec = eigs(A, nev=1, which=:SM, tol=0.0)#, maxiter=2000)
    err = sqrt(vars.hq * vars.hp) * norm(A * real.(vec) - val[1] * real.(vec))
    ψsm = reshape(real.(vec), vars.mq, vars.mp)
    normalize && ψsm ./= vars.hq * vars.hp * sum(ψsm)

    return ψsm, err
end

"Solve linear system Ax = b with reshaped output"
function linear_solver(A, b, vars::Params)
    x = A \ b
    err = norm(A * x - b)
    # x, ch = lsqr(A, b, log=true); err = ch[:resnorm][end]

    return reshape(x, vars.mq, vars.mp), err
end

function sample_gibbs(V, β, size)
    samples = zeros(size)
    for i in eachindex(samples)
        g = -π + 2π * rand()
        u = rand()
        while u > exp(-β * V(g))
            g = 2π * rand() - π
            u = rand()
        end
        samples[i] = g
    end
    return samples
end