using Parameters

function initialize_vector(range::AbstractRange)
    v = Vector{MyStruct}(undef, length(range))
    @inbounds for i in eachindex(range)
        v[i] = MyStruct(range[i], 0)
    end
    return v
end
#Come into play in multvariable case
#Dependence on num_multiparams, rstep_size, num_consts
@with_kw mutable struct MPS{R<:Real}
    rstep_size::Int64 = 1
    num_consts::Int64 = 2
    num_params::Int64 = 1# = num_multiparams - num_consts
    #m_length::Int64 = length(m)
    #a_length::Int64# = length(a)
    #alpha_length::Int64# = length(α)
    all_the_same::Bool = false# = m_length == alpha_length == a_length && alpha_length > 1
    prms_lenghts = Vector{Int64}()
    mtxtrials::Matrix{R} = zeros(R, 1, 1) #::Vector{Vector{Int64}}
    mtxparams::Matrix{R} = zeros(R, 1, 1)#::Vector{Vector{Int64}}
    num_multiparams::Int64#; @assert num_multiparams>0
end

function Base.setproperty!(mps::MPS, num_multiparams::Int64, num_consts::Int64, rand_steps::Vector{Int64},
    rstep_size::Int64 = length(rand_steps), sizes_of_parameters::Int64 ...; s::Symbol=:no_s)
    #Multiparameter separate form
     #A = [[1,2], [1,2,3], [1,4], [1], [1,5,6,7]]
    num_params::Int64 = num_multiparams - num_const
    all_the_same::Bool = false
    #To determine number of varargs
    x_tup= sizes_of_parameters
    prms_sizes::Vector{Int64} = zeros(Int64, length(x))
    for (i,u) in enumerate(a)
        j[i] =  size(u)[1]
    end
    all_the_same, p_length = all_elements_the_same(prms_sizes, check_this_size = true)
    #consts_the_same = all_elements_the_same(consts_sizes, check_this_size = true)
    multi_trialsvec::Vector{Vector{Int64}} = [[]]
    multiparamsvec::Vector{Vector{Int64}} = [[]]
    mtxtrials::Matrix{Int64}
    mtxparams::Matrix{Int64}
    for (i,k) in enumerate(rand_steps)
        multi_trialsvec[i] = k
    end
    mtxtrials = vecvec_to_array(multi_trialsvec)

    if all_the_same
        mtxparams = zeros(Int64, p_length, p_length)
    else
        for (i,k) in enumerate(rand_steps)
            multiparamsvec[i] = k
            #Convert vecvec to matrix
        end
        mtxparams = vecvec_to_array(multiparamsvec)
    end
    if s == :no_s
        print("No symbol was given to set the structure field")
    else
        error("unknown property $s")
    end
    mps(num_multiparams, rstep_size, num_const, num_params, prms_lenghts = prms_sizes, all_the_same, mtxtrials, mtxparams)
    #Let's pack all necessary information in iterable object
end
#Iterators.flatten() ; f(x) = reduce(hcat, collect(x))
#Julia 1.9 (available in beta right now) has a new stack function:
minmax(x, y) = (y < x) ? (y, x) : (x, y)
gap((min, max)) = max - min

function Base.getproperty(mps::MPS, p::Symbol)
    if p == :trials
        return getfield(mps, :mtxtrials)   # getfield() reads fields directly and it cannot be overloaded
    elseif p == :params
        return getfield(mps, :mtxparams)
    elseif p == :ncst
        return getfield(mps, :num_consts)
    elseif p == :nprs
        return getfield(mps, :num_params)
    elseif p == :rsts
        return getfield(mps, :rstep_size)
#=
    elseif p == :sls
        return getfield(mps, :SLs)   # getfield() reads fields directly and it cannot be overloaded
    elseif p == :srs
        return getfield(mps, :SRs)
    elseif p == :pls
        return getfield(mps, :PLs)
    elseif p == :prs
        return getfield(mps, :PRs)
    elseif p == :als
        return getfield(mps, :ALs)
    elseif p == :ars
        return getfield(mps, :ARs)
=#
    else
        return getfield(mps, p)
    end
end
Base.propertynames(::MPS) = (:b, :c, :q)

function set_arrays_approp(mps::MPS)
    @unpack rstep_size, num_multiparams, num_params = B(4,5,6) #get values by keys
    SLs::Array{R, 2} = zeros(R, rstep_size, NH)#This is a row as if rstep_size*NH
    SRs::Array{R, 2} = zeros(R, rstep_size, MV)
    #Multistep form
    println("num components in each vector = ", length(mps.multi_trialsvec[1]))
    println("num of vectors = ", length(mps.multi_trialsvec))

    if all_the_same
        PLLs::Array{R, 2} = zeros(R, num_params, p_length)
        PLs::Array{Float64, 2} = zeros(Float64, num_consts, p_length)
    else
        PALs::Vector{R} = zeros(R, p_length)
        #Other changable const
        #f(x) = reduce(hcat, collect(x))
        ALs::Vector{Float64} = zeros(Float64, a_length)
        ARs::Vector{Float64} = zeros(Float64, m_length)
    end
end

function fill_MPS(rand_steps, num_multiparams, num_consts; multiple_trials = false, multiple_parameters = false)
    if multiple_trials || multiple_parameters
        MPS_ = MPS()
        if multiple_parameters && length(a) == length(m)
            MPS_.pls[1] = α
            a = reshape(a, length(a), 1)#MPS_.num_consts
            m = reshape(m, length(a), 1)
            MPS_.als = hcat(a, m)
        elseif multiple_trials && rstep_size>1
            MPS_.sls[1:rstep_size, :] = rand_steps
        else
            MPS_.sls[1:rstep_size, :] = rand_steps
        end
    end
end
