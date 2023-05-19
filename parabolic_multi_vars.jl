#Symbolic representation of derivatives, etc.
#Meshes,Convex, SymEngine, CartesianGrids
using BenchmarkTools
using Base.Cartesian
using Distributions, LinearAlgebra, Distances, Calculus, SpecialFunctions, LaTeXStrings
using Printf, DataFrames, DelimitedFiles
using Symbolics, MLStyle, ToggleableAsserts
using Glob, JSON, JSONTables, CSV, Plots;
include("default_variables.jl")
include("travrse_prep.jl")
include("multips_str.jl")#p-parameters, s -steps
#include.(filter(contains(r".jl$"), readdir(pwd(); join=true)))
"""
The structure of the consistent use of the program set as follows:
1 Input data in txt files across the programm source directory
   or use default settings(Note that the flow of the programm
   depends on number of parameters&step mode defined)
2 Derived parameters are set
3 Based on number the appropriate overloads are chosen for
   traversing the grid
4 the resulting images and descriptions are set in proper
   directories for each set of parameters

    axissymmetric_parabolic_scheme(*) - solves the 1D problem in polar coordinates
    plane_parabolic_scheme(*) - solves the 2D problem in Decart coordinates

"""
#Force to garbage unused variables(in accompaniment with Var=None)
#GC.gc()
#By steps, first I am choosing variables
#read_parameters()
prpk, each_varying = init_from_files()
#Default constant(linear proportionality) and the parameter
#indicating the degree of the speed upon nonlinear dependence(Doesn't need to be so)
#m-porosity, a- constants(doesn't necessary but on the off-chance)
m, a, α = 0.1, 3, 0.3

for (k, v) in each_varying
    println("varying parameter: $k, value: $v")
end
#It's simpler to assign explicitly
m, α, a  = each_varying["m"], each_varying["α"], each_varying["a"]

all_params_number = prpk.params_num
debug_output = prpk.chosen_modes[:debug_output]
const real = prpk.chosen_modes[:real] #To avoid complex numbers in calculations
const multiple_parameters = prpk.chosen_modes[:multiple_parameters]
const multiple_trials = prpk.chosen_modes[:multiple_trials]
const debug_time_sleep =  prpk.chosen_modes[:debug_time_sleep]
formula_init = false
const save_hfig = true
same_then_reduce = true
res_dir = prpk.dir_paths[2]
fig_dir = prpk.dir_paths[3]
data_dir = prpk.dir_paths[4]
anim_dir = prpk.dir_paths[5]

used_methods = "fd"
const boundary_h = (prpk.initial_model_values["boundary_hl"], prpk.initial_model_values["boundary_hr"]) #lb::Integer= , rb::Integer=
const boundary_v = (prpk.initial_model_values["boundary_vl"], prpk.initial_model_values["boundary_vr"])
# 2D point's scope
const LD, RD, LU, RU = convert(Float64,boundary_h[1]), convert(Float64,boundary_h[2]), convert(Float64,boundary_v[1]), convert(Float64,boundary_v[2])
const hsize = RD - LD
const vsize = RU - LU
initial_shape  = prpk.initial_model_values["initial_shape"]
initial_step_value = prpk.initial_model_values["initial_step_value"]
const step_max = prpk.initial_model_values["step_max"]
const expiremntal_increase = floor(Int64, step_max/5)
_vara = 3
_varalpha = 0.3

if prpk.chosen_modes[:multiple_trials]
    Nh, Mv = 0, 0
else
    #Nh, Mv - steps without right boundary condition
    Nh, Mv = prpk.initial_model_values["Nh"], prpk.initial_model_values["Mv"]
end

rstep_size, Pdf, Pddf, mhgrid, mvgrid = choose_numgrid(hsize, vsize, initial_step_value, expiremntal_increase, step_max, Nh, Mv, a,
    multiple_trials, false)
(debug_output && println(rstep_size, Pdf, Pddf, mhgrid, mvgrid))
init_shape =" "
#From this point I am also would like to save this parameters in newly created for the form directory
#But I need number of ParametersPack's used, so I am creating beforehand all of them
#ll = mkpath(anim_dir*"$(init_shape)_ijk\\")
#print(ispath(ll))
#println( ".... $anim_dir ... $(init_shape)_ijk ... ")
#Now lets read file data or leave the default ones by my chose
initial_dir = abspath("C:/Users/2020/Desktop/CurrentReading/AsideMy/RCJsPyJulia_projects/ConcreteFluidSystem/ParabolicEquation")
if pwd() == "C:\\Users\\2020" && Sys.iswindows()
    cd(initial_dir)
    #cd("C:/Users/2020")
end
println("\nCurrent directory with program: \n", pwd())

println("Pdf:   ", eltype.(eachcol(Pdf)))
println("Pddf:   ", eltype.(eachcol(Pddf)))
#Sooner I need to write them in Result directory
CSV.write("string_primary_data.csv", Pdf)
CSV.write("model_data.csv", Pddf)
num_multiparams = prpk.params_num
mps::MPS = MPS(num_multiparams, num_consts, rand_steps,
    rstep_size = rand_steps, sizes_of_parameters)
#According to specification of Base.*, write dependence of the flow based on the overloads of functions
if multiple_parameters && multiple_trials
    write_dirs_according_to_varying_parameters(a, α, m, mhgrid, mvgrid; same_then_reduce=true, multiple_trials=true)
    #determine_initial_surface_shape(initial_shape,
elseif multiple_parameters #One-StepSet-Multi-Parameter - my case for simplicity
    write_dirs_according_to_varying_parameters(a, α, m ;same_then_reduce=true)
else
    #Must be default to One-StepSet-One-ParametersSet
    write_dirs_according_to_varying_parameters()
    #determine_initial_surface_shape(initial_shape, NL_1, NH, MV, h, τ, (LD, RD, LU, RU))
end

#CSV.write(
#------------------------------------------------------------------------------------------------------------------
#=
(debug_time_sleep && sleep(1))
Nh, Mv = Int64(Nh), Int64(Mv)
NH = Nh + 1
MV = Mv + 1



if !multiple_trials
    #Initialization of used array variables
    NH = Nh + 1 #As needed for grid including both ends
    MV = Mv + 1
    NL = zeros(2, NH)
    NL[1, 1:NH] .= 0.0
    NL[2, :] .= 0.0
    NL_1 = NL[1, :] #Passed by reference
    NL_2 = NL[2, :]
    (debug_time_sleep && println("NL: $NL, \nNL_1 and NL_2: $NL_1 and $NL_2 "))

    #Now i need to initialize shape of my surface and mesh grid also
    #SW = Point(LD, 0.0)
    #SE = Point(1.0, 0.0)
    #NW = Point(0.0, 1.0)
    #NE = Point(1.0, 1.0)
    dims = (MV, NH) # can pass in a 1-Tuple, a 2-Tuple or a 3-Tuple
    num_nodes = prod(dims)
    u = zeros(MV, NH)
    w = Nodes(Dual, dims)
    primary_data = Dict{String, Any}("used_methods" => "fd", "boundary" => Dict("boundary_h"=>Nh,"boundary_v"=>Mv),
                "number_of_steps" => Dict("NX"=>Nh, "NY"=>Mv))
    #Extract some for initialization step
    if multiple_parameters
        am = [a * j for j in a:a:num_multiparams*a]
        alpham = [α * i for i in α:α:num_multiparams*α]
        model_data = Dict("proportionality_constant" => am, "model_degree" => alpham,
                    "porosity" => m)
        MPS_.als, MPS_.pls = alpham
    else
        model_data = Dict{String, Float64}("proportionality_constant" => a, "model_degree" => α,
                    "porosity" => m)
    end
elseif !multiple_parameter
    primary_data[] = Dict("used_methods" => "fd", "boundary" => Dict("boundary_h"=>Nh,"boundary_v"=>Mv),
                "number_of_steps" => Dict("NX"=>Nh, "NY"=>Mv))
    model_data = Dict("proportionality_constant" => a, "model_degree" => α,
                "porosity" => m)
#else #Fixed step number
    #if multiple_parameter

end

#Symbolic variables used in multistep discretisized grid for
#parabolic equations
all_params_numberstr = "Params_Number_$all_params_number"
println(all_params_numberstr)
#@vars h #the unknown free surface of the fluid stream

# Generate random problem data
Symbolics.@variables A[1:Nh, 1]
A = rand(1:10, Nh, Mv)
rprimary_data, ndf = model_data_in_dir(primary_data, initial_dir)
# Create a (column vector) variable of size n x 1.
Symbolics.@variables b = randn(Mv, 1)
#grid = CartesianGrid(dims)
#To initialize h(r,t) I need the following: NH, MV, h
if !multiple_trials
    local Nh::Int64, Mv::Int64 = Pdf[3, 2]
    println("Grid size: $Nh * $Mv")
    τ, h= Pddf[1:2, 2]
    hticks, vticks= Pddf[6:7, 2]
    #Here stuck!
    NL_1= determine_initial_surface_shape!(initial_shape, multiple_trials, NL_1,
    NH, MV, h, LD, RD, LU, RU, hticks, vticks)
    initialshape_file_record = "Results/Layerh_1.txt"
    delim = " & \n"
    wvec_layer(NL_1, initialshape_file_record, 1, delim)
    plot()
end
#Shape of surface, starts with 1
delta_axis = 0.01
initial_const = 0.0

#print(u[1], "\n" , NL_1)
#println(typeof(u), "And NL_1: \n", typeof(NL_1))
(debug_time_sleep && println("U[1] type- $(u[1, :]), NL_1 type- $NL_1,\n and sizes: $(sizeof(NL_1)) ... $(sizeof(u[1, :]))"))
u[1, :]= deepcopy(NL_1)
w[1, :]= deepcopy(NL_1)
(debug_output && println("By reference?: $(NL_1 == NL[1])"))

#symmetric_layer_fill!(NL_2, NL_1, 2)
axisymmetric_parabolic_scheme(u, #timesteps::Integer,
    Pdf, Pddf, NL_1, NL_2, initial_shape)

@macroexpand @nref 2 NL i
#Additionly create plots from the resulting 2-dim massive

new_prpk = ParametersPack{Float64}()
copy!(new_prpk, prpk)#ParametersPack{Float64}(prpk.initial_model_values, params_num, control_b, pts)
println(new_prpk, "\n", prpk)
 #dst.chosen_modes
=#
