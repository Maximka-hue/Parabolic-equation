#Symbolic representation of derivatives, etc.
using SymEngine, BenchmarkTools
using CartesianGrids, Base.Cartesian #Meshes,Convex
using Distributions, LinearAlgebra, Distances, Calculus, SpecialFunctions, LaTeXStrings
using Printf, DataFrames, DelimitedFiles
using Symbolics, MLStyle, ToggleableAsserts
using Glob, JSON, JSONTables, CSV, Plots;
include("default_variables.jl")
include("travrse_prep.jl")
#include.(filter(contains(r".jl$"), readdir(pwd(); join=true)))
"""
The structure of the consistent use of the program set as follows:
1 Input data in txt files across the programm source directory
   or use default settings(Note that the flow of the programm
   depends on number of parameters&step mode defined)
2 Derived parameters are set
3 Based on number the appropriate overloads are chosen for
   traversing the grid
4 the resulting images and description are set in proper
   directories for each set of parameters

    axissymmetric_parabolic_scheme(*) - solves the 1D problem in polar coordinates
    plane_parabolic_scheme(*) - solves the 2D problem in Decart coordinates

"""
#Force to garbage unused variables(in accompaniment with Var=None)
#GC.gc()
#By steps, first I am choosing variables
#read_parameters()
prpk, each_varying = init_from_files()
for (k, v) in each_varying
    print("k: $k, v: $v")
end
new_prpk = ParametersPack{Float64}()
copy!(new_prpk, prpk)#ParametersPack{Float64}(prpk.initial_model_values, params_num, control_b, pts)
println(new_prpk, "\n", prpk)
 #dst.chosen_modes
