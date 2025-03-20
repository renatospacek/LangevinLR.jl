module Ovd2D

using Cubature
using JLLWrappers
using ForwardDiff
using Arpack
using Parameters
using LinearAlgebra
using SparseArrays
using Statistics
using IterativeSolvers

abstract type Operator end
struct FokkerPlanck <: Operator end
struct Poisson <: Operator end

include("utils_ovd2D.jl")
include("stencil_ovd2D.jl")

export Operator, FokkerPlanck, Poisson
export Params, linear_solver, spectral_solver
export get_stencil

end

