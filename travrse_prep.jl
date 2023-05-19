include("default_variables.jl")
using ToggleableAsserts, LaTeXStrings

function arithmetic_prolong(initial_step_value = 1, expiremntal_increase  = 0, step_max  = 1; av = 0, alphav = 0, num_multiparams = 1, step = true)
    if step
        rand_steps = [i for i in collect(initial_step_value:expiremntal_increase:step_max)]
    return rand_steps
    else
        #This can create the same number of parameters as grids
        am = [av * j for j in av:av:num_multiparams*av]
        alpham = [alphav * i for i in alphav:alphav:num_multiparams*alphav]
    return am, alpham
    end
end

#Compute approximate errors at tend with analytical solution
#=
function mesh_norm(u, mesh::Uniform1DFVMesh, p)
    mesh_norm(u, mesh.Δx, p)
end
function mesh_norm(u, dx::Real, p)
    @assert p > 0.0 "p must be a positive number"
    if p == Inf
        maximum(abs.(u))
    else
        (sum(abs,(u).^p)*dx)^(1/p)
    end
end
=#
#Intermittent functions
@inline function all_elements_the_same(x, check_this_size = false)
    length(x) < 2 && return true
    e1 = x[1]
    i = 2
    (check_this_size && @mayassert e1 > 1)
    @inbounds for i=2:length(x)
        x[i] == e1 || return false
    end
    return true, e1
end

function findlongest(A)
    idx = 0
    len = 0
    @inbounds for i in 1:length(A)
        l = length(A[i])
        l > len && (idx = i; len=l)
    end
    return A[idx]
end

function vecvec_to_array(vecvec)
    dim1 = length(vecvec)
    dim2 = length(vecvec[1])
    my_array = zeros(Int64, dim1, dim2)
    for i in 1:dim1
        for j in 1:dim2
            my_array[i,j] = vecvec[i][j]
        end
    end
    return my_array
end

const NArray{N} = Array{T,N} where T

function (NArray{N})(vector_of_arrays::Vector{Array{T,M}}; dims = 1) where {T, N,M}
    M <= N || throw(MethodError(NArray{N}, (vector_of_arrays,)))
    A = cat(vector_of_arrays...; dims)
    reshape(A, ntuple(i->i <= ndims(A) ? size(A,i) : 1,N))
end

"Push the derivative forward by argument h with approximation to O(h^2)"
diff_forward(f, x; l=sqrt(eps(Float64))) = (f(x+h) - f(x))/l
const DFh = diff_forward
"Advance the derivative in both directions by argument h with approximation up to O(h^3)"
diff_central(f, x; l=cbrt(eps(Float64))) = (f(x+h/2) - f(x-h/2))/l
const DCh = diff_forward
"Push the derivative backward by argument h with approximation to O(h^2)"
diff_backward(f, x; l=sqrt(eps(Float64))) = (f(x) - f(x-h))/l
const DBh = diff_forward

function spfig(initial_shape::Int64, fig_dir::String, zeros_rep::String, i::Int64, debug_output::Bool = true, save_hfig::Bool = true)
    if initial_shape in 1:3 && save_hfig
        if initial_shape == 1
            fig_path = joinpath(fig_dir, "plane_$zeros_rep$i.png")
        elseif initial_shape == 2
            fig_path = joinpath(fig_dir, "sphere_$zeros_rep$i.png")
        elseif initial_shape == 3
            fig_path = joinpath(fig_dir, "gauss_uniform_distribution_$zeros_rep$i.png")
        end
        savefig(fig_path)
    elseif !(initial_shape in 1:3)
        (debug_output && print("Released shapes are from 1 up to 3"))
    else
        plot()
    end
end


function choose_numgrid(hsize, vsize, initial_step_value, expiremntal_increase, step_max, Nh, Mv, A= 1, multiple_trials = false, multiple_parameter = false)
#Convert to positive grid steps
    hsize, vsize = abs(hsize, vsize)
    quot_stab = 2 * a[argmax(a)]
    isvy, eiy, smy = (initial_step_value^2)/quot_stab, (expiremntal_increase^2)/quot_stab, (step_max^2)/quot_stab
    rand_steps = arithmetic_prolong(initial_step_value, expiremntal_increase, step_max)
    rsy = arithmetic_prolong(isvy, eiy, smy)
    rstep_size::Int64 = 0
    if multiple_trials #Doesn't introduce new scope
        hgrid_steps::Vector{Int64} = []
        vgrid_steps::Vector{Int64} = []
    end
    if multiple_trials
         #For explicitness I check all the cases independently
        if (hsize >=1 && vsize >=1) || (hsize <= -1 && vsize  <= -1)
            hgrid_steps =  rand_steps
            vgrid_steps =  rsy
            #m- stands for multiple
            mh = hsize / hgrid_steps
            mv = vsize / vgrid_steps
        elseif hsize <= 1 && vsize <= 1 && hsize>=0 && vsize>=0
            step_max_red = convert(Int64, hsize * step_max)
            expiremntal_increase = convert(Int64, hsize * expiremntal_increase)
            rand_steps = arithmetic_prolong(initial_step_value, expiremntal_increase, step_max_red, av = 0, alphav = 0, num_multiparams = 1; step = true)
            rsy = arithmetic_prolong(isvy, expiremntal_increase^2/quot_stab, step_max_red^2/quot_stab)
            hgrid_steps = rand_steps
            vgrid_steps = rsy
            mh = hsize / hgrid_steps
            mv = vsize / vgrid_steps
        elseif hsize >= -1 && vsize >= -1 && hsize<=0 && vsize<=0
            step_max_red = convert(Int64, hsize * step_max)
            expiremntal_increase = convert(Int64, hsize * expiremntal_increase)
            rand_steps = arithmetic_prolong(initial_step_value, expiremntal_increase, step_max_red, av = 0, alphav = 0, num_multiparams = 1; step = true)
            rsy = arithmetic_prolong(isvy, expiremntal_increase^2/quot_stab, step_max_red^2/quot_stab)
            hgrid_steps = hsize* rand_steps
            vgrid_steps = vsize* rsy #Maybe twice, just check it
            mh = hsize / hgrid_steps #letter m- stands for multiple
            mv = vsize / vgrid_steps
            rstep_size = size(rand_steps)[1]
            mhs, mvs = size(mh)[1], size(mv)[1]
            @toggled_assert (length(mhs) == length(mvs) == length(rstep_size))
        else
            print("Impossible to reach out this case")
        end
        #Here data is structured into the grid parameters and specifically model ones
        Pdf = DataFrame(Names_of_parameters = ["Used Methods",
            "Boundary conditions(Cartesian)", "Number of cells along axis(By default 100*100)"],
            Params_Numbers = [used_methods, (boundary_h, boundary_v), (mhgrid, mvgrid)])
        #Here it means create them, instead using them from file
        if multiple_parameter
            #Default to the same as number of grids
            num_multiparams = rstep_size
            am, αm = arithmetic_prolong(av = a, alphav = α, num_multiparams = 1, step = false)
            @toggled_assert (length(a) == length(α) == length(rstep_size))
            Pddf = DataFrame(Names_of_derived_parameters =[
                "Time step τ", "Space step h", "Porosity m",
                "Proportionality constant a","Degree dependence α"], quantities = [mv, mh, m, am, αm])
        else
            Pddf = DataFrame(Names_of_derived_parameters =[
                "Time step τ", "Space step h", "Porosity m",
                "Proportionality constant a","Degree dependence α"], quantities = [mv, mh, m, a, α])
        end
        mhgrid = [LD:hh:RD for hh in mh]
        mvgrid = [LU:vv:RU for vv in mv]
        (hticks, vticks) = choose_approp_ticks(mhgrid, mvgrid, mh, mv, LU, LD, RU, RD)
        (debug_output && println("\n$mhgrid and $mvgrid"))
        return rstep_size, Pdf, Pddf, mhgrid, mvgrid
    else
        #All numbers of steps along also on the off-chance I am converting to abs
        h, τ = abs(hsize/(Nh+1)), abs(vsize/(Mv+1))#--------including right additional node
        hgrid = LD:h:RD # dr
        vgrid = LU:τ:RU # dr
        (hticks, vticks) = choose_approp_ticks(hgrid, vgrid, h, τ, LU, LD, RU, RD)
        (debug_output && println("\n$hgrid and $vgrid"))
        #The simplest one-set default parameters and multiprams
        Pdf = DataFrame(Names_of_parameters = ["Used Methods",
            "Boundary conditions(Cartesian)", "Number of cells along axis(By default 100*100)"],
            Params_Numbers = [used_methods, (boundary_h, boundary_v), (Nh, Mv)])
        if multiple_parameter
            num_multiparams = rstep_size
            #Types changed if a, alpha!
            am, αm = arithmetic_prolong(av = a, alphav = α, num_multiparams = 1, step = false)
            @toggled_assert (length(a) == length(α) == length(rstep_size))
            Pddf = DataFrame(Names_of_derived_parameters =[
                "Time step τ", "Space step h", "Porosity m",
                "Proportionality constant a","Degree dependence α", "hticks", "vticks"], quantities = [τ, h, m, am, αm, hticks, vticks])
        #The only set of parameters and grid configuration
        else
            Pddf = DataFrame(Names_of_derived_parameters =[
            "Time step τ", "Space step h", "Porosity m",
            "Proportionality constant a","Degree dependence α", "hticks", "vticks"],
            quantities = [τ, h, m, a, α, hticks, vticks])
        end
    end
    return 1, Pdf, Pddf, hgrid, vgrid
end

function choose_approp_ticks(hgrid::StepRange{Int64, Int64}, vgrid::StepRange{Int64, Int64},
    h::Int64, τ::Int64, LU::Int64, LD::Int64, RU::Int64, RD::Int64)
    if Nh<=30 && Mv <= 30
        hticks = hgrid
        vticks = LU:2*τ:RU
    elseif Nh<= 150 && Mv <= 100
        hticks = LD:5*h:RD
        vticks = LU:10*τ:RU
    elseif Nh<= 500  && Mv <= 500
        hticks = LD:30*h:RD
        vticks = LU:50*τ:RU
    else
        hticks = LD:40*h:RD
        vticks = LU:60*τ:RU
    end
    return hticks, vticks
end

function choose_approp_ticks(hgrid::Vector{StepRange{Int64, Int64}}, vgrid::Vector{StepRange{Int64, Int64}},
    h::Int64, τ::Int64, LU::Int64, LD::Int64, RU::Int64, RD::Int64)
    hticks::Vector{StepRange{Int64, Int64}}
    vticks::Vector{StepRange{Int64, Int64}}
    for (i, hg, vg) in enumerate(hgrid, vgrid)
        println("Checking out th $ith step!")
        hticks[i], vticks[i] = choose_approp_ticks(hg, vg,
                                                h, τ, LU, LD, RU, RD)
    end
    return hticks, vticks
end

#One-StepSet-One-ParametersSet
function determine_initial_surface_shape!(initial_shape::Int64, NL_1::Vector{Float64},
    NH::Int64, MV::Int64, h::Float64, τ::Float64, LD, RD, LU, RU,
    hticks, vticks) #Will be the ... ˜ not same, but ...
    #The number of digits to fill the name of the picture
    nd = ndigits(MV)
    zeros_rep = repeat('0', nd - 1)
    hgrid = LD:h:RD # dr
    vgrid = LU:τ:RU # dr
    if initial_shape == 1
        z::Float64 = 1.0
        @time for i in 1:NH
            NL_1[i] = z
        end
        plot(hgrid, NL_1, c= :blue, label = "h(x) = $z at timestep = 1(from first up to $Mv)", xlims= (LD, RD), xlabel="r", ylims= (LU, RU), lw=0.5,
        legend=:topleft)
        spfig(initial_shape, fig_dir, zeros_rep, 1, debug_output, save_hfig)
    elseif initial_shape == 2
        formula_init = false#use usual loops instead broadcasting
        Ru, Rd, Ld, Rd = RU^2, RD^2, LD^2, RD^2
        R = max(Ru, Rd, Ld, Rd)
        t = max(h, τ)
        if !formula_init
            initial_const = 0.0
            center = if RD- LD >= 0
                (RD + LD)/2 # RD + LD
            else
                (RD - LD)/2
            end
            NL_1[1] = initial_const
            center_step = 0
            is_center_even = false
            try
                center_step = Int(NH/2)
                is_center_even = isinteger(center_step)
            catch err
                if isa(err, InexactError)
                    println("Error with division occured")
                 #is_center_even is false
                end
            end
            if is_center_even
                center_step = Int(NH/2)
                NL_1[center_step] = R
                mid_next = center_step + 1
            else
                #center_step = round(NH/2)#the less if *.5
                center_step = floor(Int,NH/2)
                NL_1[center_step + 1] = R
                mid_next = center_step + 2
            end
            @time for i in 2:(center_step)
                print("I1: $i")
                new_value = R - ((up_step - i) * t)^2
                print("Step: $i, value: $new_value")
                #new_value = sqrt(max(0, R^2 - i * h^2))
                if new_value>=0
                    NL_1[i] = sqrt(new_value)
                else
                    NL_1[i] = 0#sqrt(-new_value)
                end
            end
            @time for j in mid_next:NH
                #print("I1: $i")
                new_value = R - ((j - up_step) * t)^2
                if new_value>=0
                    NL_1[j] = sqrt(new_value)
                else
                    NL_1[j] = 0#sqrt(-new_value)
                end
            end
            NL_1[NH] = initial_const # Or NL_1[NH] = NL[NH-1]
            plot(hgrid, NL_1, c=:blue, label=L"h(x, 0)",
            xlims= (LD, RD), xlabel="r", ylims= (LU, RU), lw=0.5,
            xticks= hticks, yticks = vticks, minorgridlinewidth = h,legend=:topleft)
            spfig(initial_shape, fig_dir, zeros_rep, 1, debug_output, save_hfig)
        else
            center_step = 0
            is_center_even = false
            try
                center_step = Int(NH/2)
                is_center_even = isinteger(center_step)
            catch err
                if isa(err, InexactError)
                    println("Error with division occured")
                 #is_center_even is false
                end
            end
            if is_center_even
               Len = Int(NH/2)
               u1 = range(LD, stop=0, length=Len)
               u2 = range(0, stop=RD, length=Len)
           else
               Len = floor(NH/2)
               u1 = range(LD, stop=0, length=Len)
               u2 = range(0, stop=RD, length=Len+1)
           end
            v = range(LU, stop=RU, length=MV)
            uu1 = -(collect(u1)^2 .- R)
            uu2 = -(collect(u2)^2 .- R)
            print(u,"\n", uu)
            NL_1 = reduce(append!, (uu1, uu2), init=Float[])
            plot(hgrid, NL_1, c=:blue, label=L"h(x, 0)",  xlims= (LD, RD), xlabel="r",
                ylims= (LU, RU), lw=0.5, xticks= hticks, yticks = vticks, minorgridlinewidth = h,
                legend=:topleft)
            spfig(initial_shape, fig_dir, zeros_rep, 1, debug_output, save_hfig)
        end
            #This defines only x axis and label refers to it's legend
            #xticks=
        elseif initial_shape == 3
            normalDensity(z) = pdf(Normal(), z)
            d0 = normalDensity.(hgrid)
            d1 = derivative.(normalDensity, hgrid)
            d2 = second_derivative.(normalDensity, hgrid)
            NL_1 = d0
            plot(hgrid, [d0 d1 d2], c=[:blue :red :green], label=[L"f(x)" L"f’(x)" L"f’’(x)"], xlims= (LD, RD), xlabel="r", ylims= (LU, RU), lw=0.5, xticks= hticks, yticks = vticks, minorgridlinewidth = h,
            legend=:topleft)
                #This defines only x axis and label refers to it's legend
            NL_1 = d0
            spfig(initial_shape, fig_dir, zeros_rep, 1, debug_output, save_hfig)
        end
        return NL_1
end

function write_dirs_according_to_varying_parameters(initial_shape::Int64, NL_1::Vector{Float64},
    NH::Int64, MV::Int64, h::Float64, τ::Float64, LD, RD, LU, RU,
    hticks, vticks, a, α, m; same_then_reduce=true)
    print("")

end

function resilient_square_root(x::Number)
        try
            sqrt(x)
        catch err
            if isa(err, DomainError)
                sqrt(complex(x))
            end
        end
    end
#=
function wvec_layer(NL_2::Vector{Float64}, file_record_p::String,
    timestep::Int64, delim = " & \n", last_time_step::Int64 = 10000, all_in_one_file = false)
    euclid_norm::Float64 = 0.0
    NS2::Int64 = 0
    if timestep>1 && timestep < last_time_step
        NS2 = size(NL_2)[1]
        euclid_norm = euclidean(NL_2, zeros(NS2))
        println("Norm of vector: $euclid_norm\t, the size of second layer: $NS2")
    end
    if isfile(file_record_p) && all_in_one_file
        # Open file in append mode and then write to it
        exampleFileIOStream =  open(file_record_p, "a")
        Base.close(exampleFileIOStream)
        #write(exampleFileIOStream, [NL_2]);
    else
        print("Opening $file_record_p")
        print(pwd())
        open(file_record_p, "w") do io
            writedlm(io, [NL_2], delim)
            if timestep>1
                writedlm(io, ["Norm of result vector: $euclid_norm",
                    "Size of NL_2: $NS2"]," & ")
            end
        end
    end
end


"""
In one parameter set I need the following:
- Two vectors for traversing the h, τ loops
- one set each of primary_data and model_data
- delimeter to write specifically results of traverses
- file path for each horizontal traverse
In multi-step configuration I would need the same, but
- additionally instead of h, τ from model_data I will
use their vector analogues
If I want to use multi-parameter case(which may be with
one-step or vector-step traverse) I will need to
choose one that will be variable throughout the cycles.
"""

"Calculations within horizontal bypass approxiamation"
function symmetric_layer_fill!(NL_2::Vector{Float64}, NL_1::Vector{Float64}, timestep::Int64)
    Rel = τ/m
    println("\n Variable steps: $τ ", " and $h"," τ/m=", "$Rel")
    println("m: $m,\tProportionality constant a: $a,\tDegree of nonlinearity: $α")
    deg_l, deg_c, deg_r = (1/α)-1, 1/α, 1+1/α
    #Boundary conditions
    #NL_2[1], maxelem_index = findmax(NL_1)
    for j in 2:Nh
         first_term = (NL_1[j+1]-NL_1[j])^deg_r/(a^deg_c* h^deg_r)
         second_term = (NL_1[j]/(a^deg_c* h))*( (NL_1[j+1]-NL_1[j])/h )^deg_c
         sech_derivative = (NL_1[j+1]-2*NL_1[j]+NL_1[j-1])/h^2
         K = sech_derivative*NL_1[j]* (NL_1[j+1]- NL_1[j])^deg_l
         (debug_output && println("Cautious at third term, : $K, others are: $second_term, $first_term"))
         third_term =  K/(α* a^deg_c* h^deg_l) #sech_derivative* NL_1[j]* (NL_1[j+1]- NL_1[j])^deg_l/(α* a^deg_c* h^deg_l)
         println("First_term: $first_term\tSecond_term: $second_term\tThird_term: $third_term\nWith second derivative = $sech_derivative")
         (debug_time_sleep && sleep(2))
         if initial_shape != 2
              NL_2[j] = Rel* (first_term- second_term- third_term) + NL_1[j]
              res_next_layer = round(NL_2[j], digits=4)
              println("NL_2[$j]: \t$res_next_layer")
         else
             res = Rel * (first_term- second_term- third_term) + NL_1[j]
             if  res < 100.0 && res > -100.0
                  NL_2[j] = res
             else
                 NL_2[j] = NL_2[(NH - j)]
            end
        end
     end
     #Non-reflective conditions
     NL_2[NH]= deepcopy(NL_2[Nh])
     NL_2[1]= deepcopy(NL_2[2])
     (debug_output && print("Second array after horizontal cross: $NL_2"))
     println("The last element approximation get used to be: $(NL_2[NH]) and first boundary condition imply at previous layer step NL_2[1] = $(NL_2[1])")
     shape_file_record = joinpath("Results", "Layerh_$(timestep).txt")
     #shape_file_record = "Results/Layerh_$(timestep).txt"
     wvec_layer(NL_2, shape_file_record, timestep)
     #norm(NL_2) norm.(eachcol(NL_2))  f1(A,d=1) = sqrt.(sum(abs2,A,d))  f1(NL_2)
     (debug_time_sleep && sleep(1))
     (debug_output && print("Last step before return: $NL_2"))
     return NL_2
end


#Let's start from the first one
function axisymmetric_parabolic_scheme(u::Array{Float64, 2}, #timesteps::Integer,
    Pdf::DataFrame, Pddf::DataFrame, uprev::Vector{Float64}, unew::Vector{Float64},
    initial_shape::Int64)
    simple_one_prmset = !multiple_trials && !multiple_parameter
    println("Pdf:   ", eltype.(eachcol(Pdf)))
    println("Pddf:   ", eltype.(eachcol(Pddf)))

    # 2D point's scope
    LD, RD, LU, RU = Float64(boundary_h[1]), Float64(boundary_h[2]), Float64(boundary_v[1]), Float64(boundary_v[2])
    #This case determine the one-parameter set with 4-point stencil
    if simple_one_prmset
        Nh::Int64, Mv::Int64 = Pdf[3, 2]
        println("Grid size: $Nh * $Mv")
        τ, h, m, a, α = Pddf[1:5, 2]#parse(Vector{Float64}
        hgrid = LD:h:RD # dr
        vgrid = LU:τ:RU # dr
        @time for i in 2:Mv
            (debug_output && print("Layer $i\t"))
            println("Before passing over to horizontal tour: layer NL_1 = ", uprev,"\nlayer NL_2 = ", unew)
            unew = symmetric_layer_fill!(unew, uprev, i)
            u[i, :] = deepcopy(unew)
            #w[i] = NL_2
            uprev = deepcopy(unew)
            unew[:] .= 0.0
            #print(hticks, vticks)
            #sleep(2)
            #, xticks= hticks, yticks = vticks, marker=:xcross, xlims= (LD, RD),ylims= (-LU, RU)
            plot(hgrid, u[i, :], c = :blue, label = "h(x) = z at timestep $i")
            plot!(hgrid, unew, label = "NL_2 at timestep $i");
                nd = ndigits(MV)
                ni = ndigits(i)
                zeros_rep = repeat('0', nd - ni)
            spfig(initial_shape, fig_dir, zeros_rep, i ,debug_output, save_hfig)
        end
    end
end
=#
