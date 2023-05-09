function arithmetic_prolong(num, av, alphav, vara = 3, varalpha = 0.3)
    if !stable_var
        am = [av * j for j in av:av:num_multiparams*av]
        alpham = [alphav * i for i in alphav:alphav:num_multiparams*alphav]
    else
        am = [vara * j for j in av:av:num_multiparams*av]
        alpham = [varalpha * i for i in alphav:alphav:num_multiparams*alphav]
    end
    return am, alpham
end

function spfig(initial_shape::Int64, fig_dir::String, zeros_rep::String, i::Int64 ,debug_output::Bool = true, save_hfig::Bool = true)
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
