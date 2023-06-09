using Parameters, Glob, LinearAlgebra, ModelingToolkit, MLStyle#Symbols in dict use @variables
#Don't forget: clear the Workspace from the Console by pressing Ctrl+D to end Julia and Enter to start it again. This works reasonably fast.
#and clearconsole() now and then
asserting() = false #making this a function results in code being invalidated and recompiled when this gets changed
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
    params_num::Int16
    #Separately from the numerical values
    control_b::Int16
    dir_paths::Vector{String}
end

#Support for inplace copying.
import Base.copy!
function copy!(dst::ParametersPack, src::ParametersPack)
    dst.initial_model_values = src.initial_model_values
    dst.chosen_modes = src.chosen_modes
    dst.params_num = src.params_num
    dst.control_b = src.control_b
    dst.dir_paths = src.dir_paths
end

ParametersPack{T}() where {T}= ParametersPack(Dict{String, T}(),  Dict{Symbol, Bool}(),
    convert(Int16, 0), convert(Int16, 0), Vector{String}())
ParametersPack{T}(iparams_num, icontrol_b, idir_paths) where {T} = ParametersPack(Dict{String, T}(),  Dict{Symbol, Bool}(),
    iparams_num, icontrol_b, idir_paths)
#ParametersPack(N, iparams_num, icontrol_b, idir_paths) = ParametersPack(Array{Float64,1}(undef, N),  Dict{Symbol, Bool}(), iparams_num, icontrol_b, idir_paths)

function variables_from_txt(txt_file::AbstractString, prpk::ParametersPack{Float64}, type_of_txt = "")
    println("Reading $txt_file")
    all_prpk = Vector{ParametersPack{Float64}}
    txt_p = uppercase(strip(type_of_txt))
    #Needed for choosing variable file from multiple
    nfilevars = 1
    each_varn = Dict{String, Vector{Float64}}()
    println("type_of_txt: $type_of_txt")
    #From this thing depends number of ParametersPack
    if txt_p== "B" || txt_p== "MODE" || txt_p== "BMODE"
        for (name, value) in zip(readdlm(txt_file, '=')[:, 1],  readdlm(txt_file, '=')[:, 2])
            name_ = lowercase(strip(name))
            println("$name_ .. $value")
            get!(prpk.chosen_modes, Symbol(name_), value)
        end
        return prpk
    #multiple_parameters define which externally file I'm using, so
    #if there already several a, m, etc. then it is already defined
    elseif txt_p== "VAR" || txt_p== "VARIABLE"
        num_of_args = Vector{Int8}()
        num_of_variable_args = Vector{Int8}() #for example m, a, α
        #num_of_args = zeros(Int8, all_parameters)
        for (n, v) in zip(readdlm(txt_file, '=')[:, 1],  readdlm(txt_file, '=')[:, 2])
            if v isa Number
                append!(num_of_args, 1)
            else
                ns, vs = size(split(n, ",") )[1], size(split(v, ",") )[1]
                println("nvars $n with(=) $ns, nvalues $v with(=) $vs")
                @mayassert ns == vs
                append!(num_of_args, ns) # parse(Int64, ns)
            end
        end
        max_param_num, ind = findmax(num_of_args)
        #Two distinct cases: values by pair and separate modes
        #max_param_num = max(readdlm(txt_file, '=')[:,2])
        for (name, value) in zip(readdlm(txt_file, '=')[:, 1],  readdlm(txt_file, '=')[:, 2])
            #Division by name
            if (size(split(name, ","))[1] == 1) && (value isa Number)
                name_ = lowercase(strip(name))
                println("Single parameter: $value, $(typeof(value))")
                get!(prpk.initial_model_values, name_, value)
            #Division by value
            elseif size(split(name, ","))[1] == 1 && size(split(value, ","))[1] > 1
                name_ = lowercase(strip(name))
                value_tmp = split(value, ",")
                #i need to parse each element in value vec
                value = parse.(Float64, value_tmp)
                print("Varying variable $name_ , $value: $(typeof(value))")
                #Just for clarity which variables I'm varying

                if name_ == "α" && (value isa Vector || value isa Tuple)
                    println("Varying α has $value ")
                    each_varn[name_] = collect(value)
                elseif name_ =="m" && (value isa Vector || value isa Tuple)
                    println("Varying m has $value ")
                    each_varn[name_] = collect(value)
                elseif name_ =="a" && (value isa Vector || value isa Tuple)
                    println("Varying a has $value ")
                    each_varn[name_] = collect(value)
                else #Multiple vars goes to the end
                    get!(prpk.initial_model_values, name_, value)
                    println(prpk)
                end
            elseif size(split(name, ","))[1] > 1 #[1] because size returns tuple
                namesm = split(name, ",")
                valuesm = split(value, ",")
                println("Multiple parameters: $namesm, $valuesm")
                for (nn, vv) in zip( namesm , valuesm )
                    vv = parse(Float64, vv)
                    #println(nn, "  -  ", vv," Expected by multiple use: ", typeof(nn), "  ", typeof(vv) )
                    get!(prpk.initial_model_values, strip(nn), vv)
                end
                   #prpk.initial_model_values[name_]
            else
                println("Parameters were mistakangly written, maybe there were zero of them")
            end
        end#for
        println("In this case, the return statement will contain a tuple")
        return (prpk, each_varn)
    elseif txt_p== "PATH" || txt_p== "PATHS"
        for  (_, value) in zip(readdlm(txt_file, '=')[:, 1],  readdlm(txt_file, '=')[:, 2])
            println(">>", value, typeof(value))
            push!(prpk.dir_paths, value);
        end
        return prpk
    else
        print("Something ... mmm, not written maybe?")
        #return nothing
    end
end

function init_from_files()
    #Default values
    result_dirsave = "C:/Users/2020/Desktop/CurrentReading/AsideMy/RCJsPyJulia_projects/ConcreteFluidSystem/ParabolicEquation"
    fig_dirsave = "C:\\Users\\2020\\Desktop\\CurrentReading\\AsideMy\\Coursework\\FlatNonNewtonFlowProblem\\Images\\IterationShapes"
    txt_parameters = read_parameters(result_dirsave, "*mode.txt")
#Reassign variables if they might be others
    paths_ = Vector{String}
    modes_ = Vector{Bool}
    variables_ = Vector{Int64}
    each_varying = Dict{String, Vector{Float64}}()
    prpk = ParametersPack{Float64}()
    #Finally I am creating coies of my ParametersPack with
    #all found varying parameters
    #for (i, (name, value)) in enumerate(each_varn)
    #    all_prpk[i] =
    for i in read_parameters(pwd(), "*modes.txt")
        if endswith(i, "modes.txt")
            println(endswith(i, "modes.txt"), "  third")
            modes_ = variables_from_txt(i, prpk, "MODE")
            println("prpk.chosen_modes  .... $prpk.chosen_modes")
        end
    end

    println("Need :multiple_parameters for the chose of variables_file, $(prpk.chosen_modes[:multiple_parameters])")
    for i in read_parameters()
       println("txt in outside loop: $i")
       if endswith(i, "pathmode.txt")
           println(endswith(i, "pathmode.txt"), "  first")
           paths_ = variables_from_txt(i, prpk, "PATH")
       elseif endswith(i, "variablesmode.txt") && prpk.chosen_modes[:multiple_parameters]
           println(endswith(i, "variablemode.txt"), "  second.m")
           variables_, each_varying = variables_from_txt(i, prpk, "VARIABLE")
           println("$variables_ ...\n $each_varying")
       elseif endswith(i, "variablemode.txt") && !prpk.chosen_modes[:multiple_parameters]
           variables_ = variables_from_txt(i, prpk, "VARIABLE")
           println("$variables_ ...\n $each_varn")
       else
           print("Reading txt_files")
       end
   end
    println("Chosen Modes will be: ", collect(values(prpk.chosen_modes)))
    control_bs, params_nums = length(prpk.chosen_modes), length(prpk.initial_model_values)
    pts = prpk.dir_paths
    #I want to set number of parameters additionally
    prpk.params_num = params_nums
    prpk.control_b = control_bs
    return prpk, each_varying
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
    #Default to all but modes file, so I additionaly indicate this point
    all_files = Vector{String}
    if read_std
    println("Standart parameters reading")
        return @match strip(pat) begin
            "" => stdrdir(dir) # only succeed when x > 5
             _ => stdrdir(dir, pat)
        end
    else
        println("Simplified reading")
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
            #myarray=int(open(readdlm,"mynums.txt"))
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

#Generate necessary input parameters for the scheme
function write_dirs_according_to_varying_parameters(a, α, m;same_then_reduce = true,  multiple_trials=false)
        init_shape = "\\"
        println(mkpath(anim_dir*"$(init_shape)_ijk\\"), ".... $anim_dir ... $(init_shape)_ijk ... ")
        for i in size(a)[1], j = size(α)[1], k = size(m)[1]
            println("Directory num: ", i, j, k)
            init_shape = @match initial_shape begin
                 1 => "Plane"
                 2 => "Sphere"
                 3 => "NormalDistribution"
                _ => "smth_sh"
            end
            if same_then_reduce && size(a)[1] == size(α)[1] == size(m)[1] && i ==j ==k
                to_create = joinpath(data_dir, "$(init_shape)_ijk")
                mkdir(to_create)
            else
                #For each a, m, α separetely
                #anim_dirn = mkpath(anim_dir*"$(init_shape)_ijk\\")
                #data_dirn = mkpath(data_dir*"$(init_shape)_ijk\\")
                to_create_anim = joinpath(anim_dir, "$(init_shape)_ijk\\")
                to_create_data = joinpath(data_dir, "$(init_shape)_ijk\\")
                (debug_output && println("Directory $to_create_anim or $to_create_data exists?"))
                println(isdir(anim_dirn) , isdir(data_dir))
                mkdir(anim_dirn)
                mkdir(data_dirn)
            end
        end
    end

#using NativeFileDialog
#path = raw"C:\Users\..."
#pick_file(joinpath.(path, "mode.txt"))
#pick_multi_file(defaultpath = path, filterlist = "txt")
#save_file(path, "txt")
