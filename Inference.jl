b"""
Parameter Inference with ABC-SMC. 
A thin wrapper around the GpABC pacakge.
"""

using DataFrames, DataStructures
using Distributions
using DocStringExtensions
using Pipe
using PyCall
using StatsBase, LsqFit

np = pyimport("numpy")
stats = pyimport("scipy.stats")
include("Utils.jl")

"""
Get mean and standard deviation for scaling purposes, accounting for `Inf`s and `0`s.
$(TYPEDSIGNATURES)
"""
function scaler(x)
    if length(x)>0
        s_mean = mean(x)
        s_std = std(x)

        if isfinite(s_mean)==false
            s_mean=0
        end

        if s_std==0
            s_std=1
        end
        return (s_mean, s_std)
    else
        return (0, 1)
    end
end

"""
Get absolute upper and lower limits for a parameter, based on the sampling distribution `d`. <br>
For uniform distributions, the limits are the lower and upper limit of the distribution. <br>
For other distributions, the 0.01st and 99.9 th percentile are calcualted. <br>
$(TYPEDSIGNATURES)
"""
function get_limits(d)
    if typeof(d)==Uniform{Float64}
        return [d.a, d.b]
    else
        lower=quantile(d, 0.001)
        upper=quantile(d, 0.999)
        return [lower, upper]
    end
end


"""
Simple Sequential ABC, loosely based on the Toni et al. algorithm.
$(TYPEDSIGNATURES)
"""
function SMC_ABC(
        priors::OrderedDict{String,Any}, 
        reference_data::Union{Array{DataFrame,1}, Array{Float64,1}, DataFrame, Float64}, 
        simulator_function::Function, 
        loss_function::Function;
        n_loops=10,
        n_samples=500,
        n_samples_init=1000,
        acceptance_quantile=0.2,
        perturbation_factor=0.05,
        progress_every=n_loops/5,
        boundaries="absorbing"
    )
    
    # infer parameter limits based on the marginal prior pdf

    limits = [get_limits(d[2]) for d in priors]
    lower = [l[1] for l in limits]
    upper = [l[2] for l in limits]
    
    # arrays to store accepted particles + distances
    accepted_particles = zeros(n_samples)
    accepted_distances = zeros(n_samples)
    
    trace_df = DataFrame()
    
    # how many samples per loop?
    n_samples_array = vcat([n_samples_init], repeat([n_samples], n_loops))
    
    # TODO: 
    # - implement cross-validation to automatically determine when to stop the calibration?
    # - also check coefficient of variation of parameter values at the end of every loop;
    #   if variation in all parameters is very low (e.g. <1%), calibration can also be stopped
    for loop in 1:(n_loops+1)
        
        if loop % progress_every==0
            @info("Loop "*string(loop)*" of "*string(n_loops+1))
        end
        
        distances = []
        particles = []
        
        # evaluate the loss function for some number of random samples
        for i in 1:n_samples_array[loop]
            let particle = zeros(length(priors))
                # in the first loop, 
                # sample from priors
                if loop==1
                    particle = [rand(x) for x in priors.vals]
                # in the 1+nth loop, 
                # sample from the particles which have been accepted in the nth loop
                else
                    particle = StatsBase.sample(accepted_particles, Weights(1 ./ accepted_distances))
                    perturbation = exp.(rand(Normal(0, perturbation_factor), size(particle)))
                    particle = particle .* perturbation
                    
                    if boundaries=="absorbing"
                        # absorbing boundaries
                        particle = max.(lower, particle)
                        particle = min.(upper, particle)
                    elseif boundaries=="none"
                        nothing
                    else
                        error("Unknwon boundary rule "*string(boundaries))
                    end
                end
                
                # generate predicted data
                predicted = simulator_function(particle)
                
                # evaluate the loss function
                let distance = Inf
                    try
                        distance = loss_function(predicted, reference_data)
                    catch
                        distance = Inf
                    end
                # store values
                append!(distances, distance)
                append!(particles, [particle])
                end
            end
        end

        # at the end of the loop:
        # filter out infinite distances/particles
        isfin = isfinite.(distances)
        
        distances = distances[isfin]
        particles = particles[isfin]
        
        if loop==1
            num_invalid=length(isfin)-sum(isfin)
            if num_invalid>1
                @info(string(num_invalid)*" of "*string(length(isfin))* " samples invalid during intialization")
            end
        end
        
        # get the best-fitting fraction of particle
        accepted_particles = particles[distances.<=quantile(distances, acceptance_quantile)]
        accepted_distances = distances[distances.<=quantile(distances, acceptance_quantile)]
        
        if length(accepted_particles)==0
            error("Did not accept any particles. Bug in loss_function or simulator_function?")
        end
        
        # record trace data
        minloc = argmin(accepted_distances)
        pt_est = accepted_particles[minloc]
        rho = accepted_distances[minloc]
        trace_df_step = DataFrame(hcat(hcat(pt_est...), rho))
        trace_df_step[!,:step] .= loop
        rename!(trace_df_step, vcat(priors.keys, "rho", "step"))
        append!(trace_df, trace_df_step)
    end
    
    accepted_df = DataFrame(vcat(accepted_particles'...))
    rename!(accepted_df, priors.keys)
    accepted_df[!,:rho] = accepted_distances
    
    return trace_df, accepted_df
end



"""
Retrieve point estimate (min distances) with marginal uncertainty (95% credible limits).
$(SIGNATURES)
"""
function get_estimate(accepted::DataFrame)
    # this is assuming that all but the last column are parameter values!
    estimate = accepted[accepted.rho.==minimum(accepted.rho),1:end-1] |> Array |> np.ravel
    CL5 = [quantile(x,.5) for x in eachcol(accepted[:,1:end-1])]
    CL95 = [quantile(x,.95) for x in eachcol(accepted[:,1:end-1])]
    return DataFrame(
        param = names(accepted)[1:end-1],
        estimate = estimate,
        CL5 = CL5,
        CL95 = CL95
    )
end

"""
Get random sample from a DataFrame of accepted particles, weighted by inverse loss.
$(SIGNATURES)
"""
function posterior_sample(accepted::DataFrame)
    weights = Weights(1 ./ accepted.rho)
    params = names(accepted)[1:end-1]
    cols = Symbol.(params)
    row_sample = sample(Array(1:nrow(accepted)), weights)
    par_sample = accepted[row_sample,cols]
    return Array(par_sample)
end

"""
Fit a dthree-parameter log logistic function to data. \n
$(SIGNATURES)
"""
function fit_log_logistic(
        df, 
        x_var,
        y_var; 
        model = @. model(x, p) = p[2]/(1+(x/p[1])^p[3])
    )
    
    # extract relevant columns
    df = @pipe df[:,[x_var, y_var]] |> 
    # drop missing
    _[completecases(_),:]
    
    try
        # set boundaries for parameters estimation
        lower_limits = [0., 0., 1e-5]
        upper_limits = [maximum(df[:,x_var]), maximum(df[:,y_var]*1.1), 20]    
        # estimate parameters + error
        fit = curve_fit(model, df[:,x_var], df[:,y_var], [mean(df[:,x_var]), mean(df[:,y_var]), 3], lower=lower_limits, upper=upper_limits)
        errors = stderror(fit)
        confidence = confidence_interval(fit, 0.05)
        # return as DataFrame
        return DataFrame(
            param = ["EC50", "upper_limit", "beta"],
            estimate = fit.param,
            std_err = errors,
            CI_05 = [max(x[1], 0) for x in confidence],
            CI_95 = [max(x[2], 0) for x in confidence]
        )
    catch
        # return missing if model can't be fitted
        return DataFrame(
            param = ["EC50", "upper_limit", "beta"],
            estimate = repeat([missing], 3),
            std_err = repeat([missing], 3),
            CI_05 = repeat([missing], 3),
            CI_95 = repeat([missing], 3)
        )
    end
end

"""
Get dose-rsponse parameters per time-point. Wrapper around ´fit_log_logistic´.
$(SIGNATURES)
"""
function dose_response_per_timepoint(
        df::DataFrame,
        x_var::Symbol, 
        y_var::Symbol)
    params = combine(groupby(df, :tday), x -> fit_log_logistic(x, x_var, y_var))
    return params
end

"""
Assign toxicity parameters. Stressor is either one of "Cu", "Ni", "Zn".
$(SIGNATURES)
"""
function assign_params(fleas_dict, particle, stressor)
    
    i = findall(x->x==stressor, ["Cu", "Ni", "Zn"])[1]
    
    fleas_dict["toxicity_particle"]["k_d"][i]    = particle[1]
    fleas_dict["toxicity_particle"]["h_max"][i]  = particle[2]
    fleas_dict["toxicity_particle"]["h_ED50"][i] = particle[3]
    fleas_dict["toxicity_particle"]["h_beta"][i] = particle[4]
    fleas_dict["toxicity_particle"]["S_max"][:,i]  = particle[5:8]
    fleas_dict["toxicity_particle"]["S_ED50"][:,i] = particle[9:12]
    fleas_dict["toxicity_particle"]["S_beta"][:,i] = particle[13:16]
    
    return fleas_dict
end

"""
Assign toxicity parameters for a particular pMoA. `i` is the index of the stressor.
$(SIGNATURES)
"""
function assign_params(fleas_dict, particle, stressor, pmoa)
    
    i = findall(x->x==stressor, ["Cu", "Ni", "Zn"])[1]
    j = findall(x->x==pmoa, ["Growth", "Maintenance", "Assimilation", "Reproduction"])[1]
    
    fleas_dict["toxicity_particle"]["S_max"] .= 0.
    
    fleas_dict["toxicity_particle"]["k_d"][i]    = particle[1]
    fleas_dict["toxicity_particle"]["h_max"][i]  = particle[2]
    fleas_dict["toxicity_particle"]["h_ED50"][i] = particle[3]
    fleas_dict["toxicity_particle"]["h_beta"][i] = particle[4]
    fleas_dict["toxicity_particle"]["S_max"][j,i]  = particle[5]
    fleas_dict["toxicity_particle"]["S_ED50"][j,i] = particle[6]
    fleas_dict["toxicity_particle"]["S_beta"][j,i] = particle[7]
    
    return fleas_dict
end