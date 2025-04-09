module Lang1D

using Cubature
using JLLWrappers
using ForwardDiff
using Arpack
using Parameters
using LinearAlgebra
using SparseArrays
using Statistics
using IterativeSolvers
using QuadGK

using ..Common: Operator, FokkerPlanck, Poisson, meshgrid
export Operator, FokkerPlanck, Poisson, meshgrid

include("utils_lang1D.jl")
include("stencil_lang1D.jl")

export Params, linear_solver, spectral_solver
export get_stencil

end

