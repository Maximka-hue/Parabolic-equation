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

"Push the derivative forward by argument h with approximation to O(h^2)"
diff_forward(f, x; l=sqrt(eps(Float64))) = (f(x+h) - f(x))/l
const DFh = diff_forward
"Advance the derivative in both directions by argument h with approximation up to O(h^3)"
diff_central(f, x; l=cbrt(eps(Float64))) = (f(x+h/2) - f(x-h/2))/l
const DCh = diff_forward
"Push the derivative backward by argument h with approximation to O(h^2)"
diff_backward(f, x; l=sqrt(eps(Float64))) = (f(x) - f(x-h))/l
const DBh = diff_forward


function choose_numgrid(hsize, vsize, initial_step_value, expiremntal_increase, step_max, Nh, Mv, multiple_trials = false, multiple_parameter = false)
    rstep_size::Int64 = 0
    if multiple_trials #Doesn't introduce new scope
        hgrid_steps::Vector{Int64} = []
        vgrid_steps::Vector{Int64} = []
    end
    if multiple_trials
        if hsize >=1 && vsize >=1
            rand_steps = arithmetic_prolong(initial_step_value, expiremntal_increase, step_max)
            hgrid_steps = hsize* rand_steps
            vgrid_steps = vsize* rand_steps
            #m- stands for multiple
            mh = hsize / hgrid_steps
            mv = vsize / vgrid_steps
        elseif hsize <= 1 && vsize <= 1 && hsize>=0 && vsize>=0
            step_max_red = convert(Int64, hsize* step_max)
            expiremntal_increase = step_max_red
            rand_steps = arithmetic_prolong(initial_step_value, expiremntal_increase, step_max,av = 0, alphav = 0, num_multiparams = 1; step = true)
            rand_steps = arithmetic_prolong(initial_step_value, expiremntal_increase, step_max)
            hgrid_steps = rand_steps
            vgrid_steps = rand_step
            mh = hsize / hgrid_steps
            mv = vsize / vgrid_steps
        elseif hsize <= -1 && vsize  <= -1
            #Convert to positive grid steps
            hsize = min(0,- hsize)
            vsize = min(0,- vsize)
            rand_steps = arithmetic_prolong(initial_step_value, expiremntal_increase, step_max)
            hgrid_steps = hsize* rand_steps
            vgrid_steps = vsize* rand_steps
            mh = hsize / hgrid_steps
            mv = vsize / vgrid_steps
        elseif hsize >= -1 && vsize >= -1 && hsize<=0 && vsize<=0 #For explicitness
            hsize = min(0,- hsize)
            vsize = min(0,- vsize)
            step_max_red = hsize* step_max
            expiremntal_increase = step_max_red
            rand_steps = arithmetic_prolong(initial_step_value, expiremntal_increase, step_max)
            hgrid_steps = rand_steps
            vgrid_steps = rand_step
            mh = hsize / hgrid_steps #m- stands for multiple
            mv = vsize / vgrid_steps
            rstep_size = size(rand_steps)[1]
            mhs, mvs = size(mh)[1], size(mv)[1]
            @toggled_assert (length(mhs) == length(mvs) == length(rstep_size))
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
        (debug_output && println("\n$mhgrid and $mvgrid"))
        return rstep_size, Pdf, Pddf, mhgrid, mvgrid
    else
        if hsize>=0 && vsize>=0
            h = (hsize/Nh)
            τ = (vsize/Mv)
        else
            h = -(hsize/Nh)
            τ = -(vsize/Mv)
        end
        hgrid = LD:h:RD # dr
        vgrid = LU:τ:RU # dr
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
                "Proportionality constant a","Degree dependence α"], quantities = [τ, h, m, am, αm])
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
