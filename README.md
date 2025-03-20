# LangevinLR.jl

A Julia implementation of some PDE-based methods for computing the linear response (and more generally, nonlinear/full response curves) for underdamped and overdamped Langevin dynamics.

The approaches here implemented rely on finite-difference discretizations of the relevant differential operator (e.g., generator of the associated dynamics), used to numerically approximate the solution to the PDE at hand.

## Approaches

### Mathematical description

There are currently two PDE approaches implemented for computing the (linear) response. Letting $\mathcal{L}$ denote the generator of the associated stochastic dynamics, and $\mathcal{L}^\dagger$ its $L^2(\mathcal{X})$-adjoint, the methods are as follows:

1. **Poisson equation (equilibrium method):** The linear response $\rho$ of some observable $f$ can be computed as 

$$ \rho = \langle \phi, f\rangle_{L^2(\mu)} = \int_\mathcal{X} \phi f \, d\mu,$$

with $\phi$ the solution to the Poisson equation $-\mathcal{L}\phi = f$.

Due to being an equilibrium-based approach, this method is restricted to linear responses only.

2. **Fokker-Planck (nonequilibrium method):** Alternatively, $\rho$ can be computed by solving for the nonequilibrium invariant measure $\mu_\eta$ directly:

```
```math 
\rho = \lim_{\eta\to 0} \frac{1}{\eta}\mathbb{E}_\eta(R) = \lim_{\eta\to 0} \frac{1}{\eta}\int_\mathcal{X} f \, d\mu_\eta,
```
```

with $\mu_\eta$ the solution to the stationary Fokker-Planck equation $\mathcal{L}^\dagger\mu_\eta = 0$. More generally, this method allows for computing the full nonlinear response profile $\mathbb{E}_\eta(R)$ for arbitrary values of $\eta$.

### Numerical approximation

Both methods can be numerically realized by discretizing the PDE operator at hand to solve the linear system. The (linear) response can then be obtained with a simple quadrature.

See `example_ovd2D.jl` for an example of both implementations with a detailed step-by-step explanation.

## Installation

To install `LangevinLR`, add it directly from source in the Pkg REPL:

```bash
] add https://github.com/renatospacek/LangevinLR.jl
```
