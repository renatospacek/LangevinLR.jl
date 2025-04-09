using Plots
using Printf
using LangevinLR.Lang1D

V(x) = (1.0 - cos(x)) / 2.0

force(y) = 1.0
force() = force(0.0)

# ==========================================================================================
function main()
    mq = 51
    mp = 200
    Lp = 5.0
    vars = Params(Lp, mq, mp, 1.0, 1.0, V, force)
    η = 0.1

    f = force() * vars.Y

    L_FP = get_stencil(vars, FokkerPlanck(), η)
    ψ_η, err_FP = spectral_solver(L_FP, vars)

    lr_FP = vars.hq * vars.hp * sum(f .* ψ_η) / η

    L_poisson = get_stencil(vars, Poisson())
    ϕ, err_poisson = linear_solver(-L_poisson, vcat(f...), vars)

    ψ_0 = vars.gibbs.(vars.X, vars.Y)
    lr_poisson = vars.hq * vars.hp * sum(ϕ .* f .* ψ_0)

    @printf("Poisson residual = %.4g \n", err_poisson)
    @printf("Fokker-Planck residual = %.4g \n", err_FP)
    @printf("Poisson linear response = %.4g \n", lr_poisson)
    @printf("Fokker-Planck linear response = %.4g \n", lr_FP)

    ψ_0 = vars.gibbs.(vars.X, vars.Y) # discrete Gibbs distribution of size m by m

    plt_FP = surface(vars.X, vars.Y, ψ_η, legend=nothing, title="ψ_η (η = $η)", xlabel="q", ylabel="p")
    plt_gibbs = surface(vars.X, vars.Y, ψ_0, legend=nothing, title="Gibbs dist. ψ_0", xlabel="q", ylabel="p")
    plt = plot(plt_gibbs, plt_FP, layout=2)
    display(plt)
end

@time main()
