using Plots
using Printf
using LangevinLR.Ovd2D
#=
    ##### Outline of how to compute the linear response, step-by-step:

    1) Define your potential `V` below
        (if you're doing nonequilibrium, also define `force`)
    2) Instantiate `Params`, which contains simulation parameters
    3) Discretize operator. There are 2 options here:
        3a) Solve Poisson equation -Lϕ = f, with L the generator of the dynamics
            (so linear response = <ϕ,f> = ∫_X ϕf ψ_0, with ψ_0 the gibbs distribution)[^1]
        3b) Solve Fokker-Planck L^† ψ_η = 0 for the nonequilibrium invariant measure ψ_\eta,
                with L^† the L^2(Ω)-adjoint of L
            (so linear response = \lim_{η→0} 1/η \int_Ω f ψ_η)
        
        The discretization of L and L^† can be obtained with `get_stencil()` as follows:
        
            L = get_stencil(vars, Poisson())
            L^† = get_stencil(vars, FokkerPlanck())

    4) Solve the corresponding linear system:
        4a) For Poisson, use `linear_solver(L, f)`
            L is an m^2 by m^2 matrix, and f is an m^2 vector. You can use `vcat(x...)` to 
            flatten a matrix into a vector (see examples below)
        4b) For FP, use `spectral_solver(L^†, vars)`
    
    [^1] For overdamped it's actually 1 - <ϕ,f> (see 10.1093/imanum/dru056 proof of Lemma 3.2)
=#

###----- Define your potential (and force) here -----###
V(x, y) = cos(2π*x) + cos(4π*y) + 0.5*cos(2π*x)*cos(4π*y)

force(y) = [1.0, 0.0]
force() = force(0.0)

# ==========================================================================================
function main()
    ##### Parameters
    #=
        - `vars` holds general simulation parameters. Initialize it with:
            `L`: length of the domain [0, L] (depends on choice of potential `V`)
            `m`: number of discretization points per dimension
            `β`: inverse temperature
            `V`: potential energy function
            `force`: perturbation forcing
        - η only needs to be defined if you're doing nonequilibrium
            (i.e., it's set to 0 by default if not defined)
    =#
    L = 1.0
    m = 100
    β = 1.0
    η = 0.1
    vars = Params(L, m, β, V, force)
    
    #=
        - in this example, we compute the mobility/diffusivity, for which the observable of 
        interest is ∇V in the direction of the perturbation, i.e., Fᵀ∇V
        - this observable corresponds to `f`, i.e., the right-hand side of the Poisson eq as well as 
        the quantity integrated wrt ψ_η for the nonequilibrium approach
    =#
    obs(x,y) = force()[1].*vars.d₁V(x,y) + force()[2].*vars.d₂V(x,y)
    f = obs.(vars.Q1, vars.Q2)

    ##### Example 1: Solving Poisson equation -Lϕ = f --------------------------------------
    #=
        - get discretized generator `L_poisson` from `get_stencil(vars, Poisson(), η)`
            (only supports equilibrium, i.e., η = 0 [ommitting η defaults to 0])
        - solve linear system Ax = b with `x = linear_solver(A, b)`
            (note that linear_solver reshapes `x` and returns it as an m by m matrix)
    =#
    L_poisson = get_stencil(vars, Poisson())
    ϕ, err_poisson = linear_solver(L_poisson, vcat(f...))

    ### Computing the mobility with Poisson solution
    # linear response given by <ϕ,f> = ∫_X ϕf ψ_0 
    ψ_0 = vars.gibbs.(vars.Q1, vars.Q2) # discrete Gibbs distribution of size m by m
    lr_poisson = vars.h^2*sum(ϕ.*f.*ψ_0) # quadrature yields linear response

    ##### Example 2: Solving Fokker-Planck equation L^† ψ_η = 0 ----------------------------
    #=
        - get discretized generator `L^†` from `get_stencil(vars, FokkerPlanck(), η)`
      ``      (if η not passed, it's treated as 0, so make sure to specify η here)
        - `x = spectral_solver(A, vars)` solves for x in Ax = 0
    =#
    L_FP = get_stencil(vars, FokkerPlanck(), η)
    ψ_η, err_FP = spectral_solver(L_FP, vars)

    ### Computing the mobility with Fokker-Planck
    # linear response given by 1/η \int_X f ψ_η
    lr_FP = vars.h^2*sum(f.*ψ_η)/η 

    @printf("Poisson residual = %.4g \n", err_poisson)
    @printf("Fokker-Planck residual = %.4g \n", err_FP)
    @printf("Poisson linear response = %.4g \n", lr_poisson)
    @printf("Fokker-Planck linear response = %.4g \n", lr_FP)

    ### visualize ψ_η and ψ_0
    # plt_FP = surface(vars.Q1, vars.Q2, ψ_η, legend=nothing, title="ψ_η (η = $η)")
    # plt_gibbs = surface(vars.Q1, vars.Q2, ψ_0, legend=nothing, title="Gibbs dist. ψ_0")
    # plt = plot(plt_gibbs, plt_FP, layout = 2)
    # display(plt)
end

@time main()
