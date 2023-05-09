using Parameters, Glob, LinearAlgebra, ModelingToolkit#Symbols in dict use @variables
#Don't forget: clear the Workspace from the Console by pressing Ctrl+D to end Julia and Enter to start it again. This works reasonably fast.
#and clearconsole() now and then
asserting() = true #making this a function results in code being invalidated and recompiled when this gets changed
const all_parameters = 12

macro mayassert(test)
  esc(:(if $(@__MODULE__).asserting()
    @assert($test)
   end))
end


mutable struct ParametersPack{T}
    initial_model_values::Dict{String,T}
    #With this dictionary from my point of view will be
    #better to collect several possible strings
    #into the unique symbol
    chosen_modes::Dict{Symbol, Bool}
    params_num::UInt16
    #Separately from the numerical values
    control_b::UInt16
    dir_paths::Vector{String}
end

ParametersPack{Float64}() where {T}= ParametersPack(Dict{String, Float64}(),  Dict{Symbol, Bool}(), 0, 0, Vector{String}())
ParametersPack{T}(iparams_num, icontrol_b, idir_paths) = ParametersPack(Dict{String, T}(),  Dict{Symbol, Bool}(), iparams_num, icontrol_b, idir_paths)
#ParametersPack(N, iparams_num, icontrol_b, idir_paths) = ParametersPack(Array{Float64,1}(undef, N),  Dict{Symbol, Bool}(), iparams_num, icontrol_b, idir_paths)

function variables_from_txt(txt_file::AbstractString, prpk::ParametersPack{Float64}, type_of_txt = "")
    println("Reading $txt_file")
    all_prpk = Vector{ParametersPack{Float64}}
    txt_p = uppercase(strip(type_of_txt))
    #Needed for choosing variable file from multiple
    nfilevars = 1
    each_varn = Dict{String, Vector{Float64}}()
    #From this thing depends number of ParametersPack
    if txt_p== "B" || txt_p== "MODE" || txt_p== "BMODE"
        for name, value in readdlm(txt_p, '=')[:,1], readdlm(txt_p, '=')[:,2]
            name_ = lowercase(strip(name))
            get!(prpk.chosen_modes, name_, value)
            return prpk
    end
    #multiple_parameters define which externally file I'm using, so
    #if there already several a, m, etc. then it is already defined
    elseif txt_p== "VAR" || txt_p== "VARIABLE"
        num_of_args = Vector{Int8}()
        num_of_variable_args = Vector{Int8}() #for example m, a, α
        #num_of_args = zeros(Int8, all_parameters)
        for (n, v) in (readdlm(txt_p, '=')[1,:], readdlm(txt_p, '=')[2,:])
            if v isa Number
                append!(num_of_args, 1)
            else
                ns, vs = size(split(n, ","))[1], size(split(v, ","))[1]
                println("nvars $n with(=) $ns, nvalues $v with(=) $vs")
                @mayassert ns == vs
                append!(num_of_args, ns)
            end
        end
        max_param_num, ind = findmax(num_of_args)
        #Two distinct cases: values by pair and separate modes
        max_param_num = max(readdlm(txt_p, '=')[:,1]
        for name, value in readdlm(txt_p, '=')[:,1], readdlm(txt_p, '=')[:,2]
            if size(split(name, ","))[1] > 1 #[1] because size returns tuple
                println("Multiple parameters")
                for nn, vv in split(name, ","), split(value, ",")
                    get!(prpk.initial_model_values, nn, vv)
                    #Just for clarity which variables I'm varying
                    if name_ =="α" && (value isa Vector || value isa Tuple)
                        println("Varying α has $value ")
                        each_varn[name_] = value
                    elseif name_ =="m" && (value isa Vector || value isa Tuple)
                        println("Varying m has $value ")
                        each_varn[name_] = value
                    elseif name_ =="a" && (value isa Vector || value isa Tuple)
                        println("Varying a has $value ")
                        each_varn[name_] = value
                    end
                end
                   #prpk.initial_model_values[name_]
            elseif size(split(name, ","))[1] = 1
                println("Single parameter")
                get!(prpk.initial_model_values, name, value)
            else
                println("Parameters were mistakangly written, maybe there were zero of them")
            end
        end
        return (prpk, each_varn)
    elseif txt_p== "PATH" || txt_p== "PATHS"
        for _, value in readdlm(txt_p, '=')[:,1], readdlm(txt_p, '=')[:,2]
            prpk.dir_paths.push!(value)
        end
        return prpk
    else
        print("Something ... mmm, not written maybe?")
        return nothing
    end
end
#Support for inplace copying.
import Base.copy!
function copy!(dst::ParametersPack, src::ParametersPack)
    dst.initial_model_values .= src.initial_model_values
    dst.chosen_modes = src.chosen_modes
end

function copy!(dst::ParametersPack, src::Dict{String, Bool})
    #Default values
    result_dirsave = "C:/Users/2020/Desktop/CurrentReading/AsideMy/RCJsPyJulia_projects/ConcreteFluidSystem/ParabolicEquation"
    fig_dirsave = "C:\\Users\\2020\\Desktop\\CurrentReading\\AsideMy\\Coursework\\FlatNonNewtonFlowProblem\\Images\\IterationShapes"
    txt_parameters = read_parameters(result_dirsave, "*mode.txt")
#Reassign variables if they might be others
    paths_ = Vector{String}
    modes_ = Vector{Bool}
    variables_ = Vector{Int64}
    prpk = ParametersPack(N, iparams_num, icontrol_b, idir_paths)
    #Finally I am creating coies of my ParametersPack with
    #all found varying parameters
    #for (i, (name, value)) in enumerate(each_varn)
    #    all_prpk[i] =
    for i in read_parameters()
       println("$i")
       if endswith(i, "pathmode.txt")
           paths_ = variables_from_txt(i, prpk, type_of_txt = "PATH")
       elseif endswith(i, "variablemode.txt")
           variables_, each_varn = variables_from_txt(i, prpk, type_of_txt = "VARIABLE")
           println("$variables_ ... $each_varn")
       elseif endswith(i, "mode.txt")
           modes_ = variables_from_txt(i, prpk, type_of_txt = "MODE")
       else
           print("Reading txt_files")
       end
   end
    ml_prms,#multiple_parameters about parameters like α, m, a
        deb_ts,#additional sleep time beetween function evaluations
        deb_prout,#debug_output additional output
        real,#type stable evaluation witout complex numbers
        ml_trials,#define step modified repetition
        save_hfig#(h- horizontal)the name speaks for itself
        = @match read_parameters() begin

    ParametersPack(N, iparams_num, icontrol_b, Vector[result_dirsave, fig_dirsave])
    dst.chosen_modes
end


#------------------------------------------------------------------------------
#Simply for scaning data files
function rdir(dir::AbstractString, pat::Glob.FilenameMatch = fn"*.txt")
    result = String[]
    for (root, dirs, files) in walkdir(dir)
        append!(result, filter(contains(pat), readdir(dir; join=true)),#filter!(f -> occursin(pat, f), files),#all matching pattern
            joinpath.(root, files))
    end
    return result#i.e. all files like String
end
#readdir(dir::AbstractString=pwd();
      #join::Bool = false, # <-- HERE
      #sort::Bool = true,
  #) -> Vector{String})    -----------------function is in Base

function stdrdir(dir::AbstractString, pat::Glob.FilenameMatch = fn"*.txt") # = fn".txt$")
    result = filter(contains(pat), readdir(dir; join=true))#all matching pattern
    return result#i.e. all files like String
end

rdir(dir::AbstractString, pat::Glob.FilenameMatch) = rdir(dir, Glob.FilenameMatch(pat))
stdrdir(dir::AbstractString, pat::AbstractString) = stdrdir(dir, Glob.FilenameMatch(pat)) #the same as above

function read_parameters(dir::AbstractString  = pwd(), pat::AbstractString = raw"*mode.txt", read_std = true)
    all_files = Vector{String}
    if read_std
    print("Standart parameters reading")
        return @match strip(pat) begin
            "" => stdrdir(dir) # only succeed when x > 5
             _ => stdrdir(dir, pat)
        end
    else
        print("Simplified reading")
        return @match strip(pat) begin
            "" => @time rdir(dir) # only succeed when x > 5
             _ => @time rdir(dir, pat)
        end
    end
end

#Another approach
function model_data_in_dir(primary_data, initial_dir)
    # pass data as a json string (how it shall be displayed in a file)
    string_primary_data = JSON.json(primary_data)
    # write the file with the stringdata variable information
    all_data_files = rdir(initial_dir, "*data*.txt")
    #sort(readdir("."); by=f->stat(f).inode)
    #df = glob("data*.txt"i) || glob("*data*.txt"i)
    rprimary_data = Dict()
    ndf = size(all_data_files)[1]
    println("The acquired files: $all_data_files, \nWith size: $ndf")
    if ndf == 1
        fn = all_data_files[1]
        println("File found in the current directory: $fn")
        open(fn, "r") do f
            global rprimary_data
            dicttxt = readall(all_data_files)  # file information to string
            rprimary_data= JSON.parse(dicttxt)  # parse and transform data
        end
    elseif ndf == 0
        println("There are no files in directory")
    elseif ndf >=1
        for (i, file) in enumerate(all_data_files)
             mv(file, "$file_$i.txt")
        end
    end
    return (rprimary_data, ndf)
end
#ParametersPack(mode::Dict{Symbol, Bool}) =


#Generate necessary input parameters for the scheme


used_methods::AbstractString = "fd"
const boundary_h = (-5, 5) #lb::Integer= , rb::Integer=
const boundary_v = (-5, 5)
hgrid_steps::Vector{Int64} = []
vgrid_steps::Vector{Int64} = []
const Nh, Mv::Int64= 50, 100
initial_shape = 3
formula_init = false
# 2D point's scope
const LD, RD, LU, RU = Float64(boundary_h[1]), Float64(boundary_h[2]), Float64(boundary_v[1]), Float64(boundary_v[2])
const hsize = RD - LD
const vsize = RU - LU
const initial_step_value::Int64 = 10
const step_max::Int64 = 1000
const expiremntal_increase = floor(Int64, step_max/5)
#constant(linear proportionality) and the parameter
#indicating the degree of the speed upon nonlinear dependence
a, α = 3, 1.2
#Convinience parameter and specifying porosity(constant (porosity))
m = 1

#using NativeFileDialog
#path = raw"C:\Users\..."
#pick_file(joinpath.(path, "mode.txt"))
#pick_multi_file(defaultpath = path, filterlist = "txt")
#save_file(path, "txt")
