
const real = true
const multiple_parameter = false
const debug_time_sleep = false
const debug_output = true
const multiple_trials = false
const save_hfig = true
fig_dir = "C:\\Users\\2020\\Desktop\\CurrentReading\\AsideMy\\Coursework\\FlatNonNewtonFlowProblem\\Images\\IterationShapes"

#Generate necessary input parameters for the scheme

const Params_Number = 5
used_methods::AbstractString = "fd"
const boundary_h = (-5, 5) #lb::Integer= , rb::Integer=
const boundary_v = (-5, 5)
hgrid_steps::Vector{Int64} = []
vgrid_steps::Vector{Int64} = []
const Nh, Mv::Int64= 50, 50
# 2D point's scope
const LD, RD, LU, RU = Float64(boundary_h[1]), Float64(boundary_h[2]), Float64(boundary_v[1]), Float64(boundary_v[2])
const hsize = RD - LD
const vsize = RU - LU
const initial_step_value::Int64 = 10
const step_max::Int64 = 1000
const expiremntal_increase = floor(Int64, step_max/5)
