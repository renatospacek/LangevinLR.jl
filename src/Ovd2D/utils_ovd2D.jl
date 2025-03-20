# TODO generalize domains; currently only of the form [0,L)
"Construct 2D meshgrid"
function meshgrid(x, y)
    X = x' .* ones(length(y))
    Y = ones(length(x))' .* y

    return X, Y
end

"""
Struct for general parameters

    - `L`: length of the q domain: q ∈ [0, L]
        (needs to be adjusted if potential is changed)
    - `m`: number of discretization points of q
    - `h`: spatial stepsize
    - `q`: discretized vector
    - `Q1, Q2`: m by m grids of q

    NOTE: current implementation assumes the periodic domain is [0, L). For domains centered
    around 0, simply shift `q` appropriately
"""
@with_kw mutable struct Params
    L::Float64
    m::Int64
    h::Float64
    q::Array{Float64,1}
    Q1::Array{Float64,2}
    Q2::Array{Float64,2}
    β::Float64
    V::Function
    d₁V::Function
    d₂V::Function
    ΔV::Function
    force::Function
    gibbs::Function

    function Params(L, m, β, V, force)
        h = L/m
        q = range(0.0, L - h + 1e-12, m)
        Q1, Q2 = meshgrid(q, q)

        d₁V, d₂V, ΔV = diff_potential(V)
        gibbs = _get_gibbs(V, L, β)
        
        new(L, m, h, q, Q1', Q2', β, V, d₁V, d₂V, ΔV, force, gibbs)
    end
end

"Compute partial derivatives and Laplacian for a given potential `V`"
function diff_potential(V::Function)
    d₁V(x, y) = ForwardDiff.derivative(x -> V(x, y), x)
    d₂V(x, y) = ForwardDiff.derivative(y -> V(x, y), y)
    dd₁V(x, y) = ForwardDiff.derivative(x -> d₁V(x, y), x)
    dd₂V(x, y) = ForwardDiff.derivative(y -> d₂V(x, y), y)
    ΔV(x, y) = dd₁V(x, y) + dd₂V(x, y)

    return d₁V, d₂V, ΔV
end

"Define method for Gibbs distribution"
function _get_gibbs(V::Function, L::Float64, β::Float64)
    Z = hcubature(x -> exp(-β*V(x...)), [0, 0], [L, L])[1]

    return ((x, y) -> exp(-β*V(x, y))/Z)
end

"Solve Ax = 0 for sparse A by finding eigenvec associated with the 0 eigenvalue"
function spectral_solver(A, vars::Params; normalize=true)
    # @unpack_Params vars

    val, vec = eigs(A, nev=1, which=:SM, tol=0.0)#, maxiter=2000)
    err = sqrt(vars.h)*norm(A*real.(vec) - val[1]*real.(vec))
    ψsm = reshape(real.(vec), vars.m, vars.m)
    normalize && ψsm ./= vars.h^2*sum(ψsm)

    return ψsm, err
end

"Solve linear system Ax = b with reshaped output"
function linear_solver(A, b)
    m = Int(sqrt(length(b)))

    x = A\b; err = norm(A*x - b)
    # x, ch = lsqr(A, b, log=true); err = ch[:resnorm][end]

    return reshape(x, m, m), err
end
