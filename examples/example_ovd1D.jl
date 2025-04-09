using Revise
using Plots
using Printf
using LangevinLR.Ovd1D
plotlyjs()

V(x) = (1.0 - cos(x)) / 2.0

force(y) = 1.0
force() = force(0.0)

# ==========================================================================================
function main()
    L = pi
    m = 100
    β = 1.0
    vars = Params(L, m, β, V, force)
    η = 1.5

    f = force() * vars.dV.(vars.q)

    L_FP = get_stencil(vars, FokkerPlanck(), η)
    ψ_η, err_FP = spectral_solver(L_FP, vars)

    lr_FP = vars.h * sum(f .* ψ_η) / η

    L_poisson = get_stencil(vars, Poisson())
    ϕ, err_poisson = linear_solver(L_poisson, f)
    ψ_0 = vars.gibbs.(vars.q)

    lr_poisson = vars.h * sum(ϕ .* f .* ψ_0)

    @printf("Poisson residual = %.4g \n", err_poisson)
    @printf("Fokker-Planck residual = %.4g \n", err_FP)
    @printf("Poisson linear response = %.4g \n", lr_poisson)
    @printf("Fokker-Planck linear response = %.4g \n", lr_FP)

    ψ_0 = vars.gibbs.(vars.q) # discrete Gibbs distribution of size m by m

    plt_FP = plot(vars.q, ψ_η, legend=nothing, title="ψ_η (η = $η)")
    plt_gibbs = plot(vars.q, ψ_0, legend=nothing, title="Gibbs dist. ψ_0")
    plt = plot(plt_gibbs, plt_FP, plt_diff, layout=2)
    display(plt)
end

@time main()
