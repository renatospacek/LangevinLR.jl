module Ovd1D

using QuadGK
using JLLWrappers
using ForwardDiff
using Arpack
using Parameters
using LinearAlgebra
using SparseArrays
using Statistics
using IterativeSolvers

using ..Common: Operator, FokkerPlanck, Poisson, meshgrid
export Operator, FokkerPlanck, Poisson, meshgrid

include("utils_ovd1D.jl")
include("stencil_ovd1D.jl")

export Params, linear_solver, spectral_solver
export get_stencil

end

