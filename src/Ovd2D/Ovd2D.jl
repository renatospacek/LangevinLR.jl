module Ovd2D

using StaticArrays
using Cubature
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

include("utils_ovd2D.jl")
include("stencil_ovd2D.jl")

export Params, linear_solver, spectral_solver
export get_stencil

end

