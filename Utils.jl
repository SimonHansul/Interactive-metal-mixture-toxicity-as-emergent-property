# Convenience functions
using PyCall
using DocStringExtensions
using Pipe
using CSV
py"""
import numpy as np
"""

add_dim(x::Array) = reshape(x, (size(x)...,1))

"""
drop missing values in pipe
$(SIGNATURES)
"""
function drop_na(df::DataFrame)
    return df[completecases(df),:]
end

function drop_na(ar::Array{Union{Missing,Float64},1})
    mask1 = .!ismissing.(ar)
    mask2 = isfinite.(ar)
    return ar[mask1.&mask2]
end

"""
Enumerate unique values
"""
function enique(x)
    return enumerate(unique(x))
end

"""
Read observed data
$(SIGNATURES)
"""
function read_observed(path::String; timevar=:tday, yvar=:count)
    observed = @pipe CSV.read(path, DataFrame) |> leftjoin(_, combine(groupby(_[_.treatment.=="Co",:], timevar), yvar => mean), on=timevar)
    observed.response = observed[:,yvar] ./ observed[:,Symbol(string(yvar)*"_mean")]
    observed.tday = convert.(Float64, observed.tday)
    return observed
end


"""
Check if timepoint in timepoints
"""
function t_in_timepoints(x)
    return x in timepoints
end

"""
Plot population trajectories in different treatments, given predicted and observed data.
$(SIGNATURES)
"""
function plot_all(predicted::DataFrame, observed::DataFrame; alpha=.2, co_mean=nothing)
    trts = unique(predicted.treatment)
    fig, ax = plt.subplots(ncols=size(trts)[1], figsize=(20, 4), sharey=true)
    for (i, trt) in enumerate(trts)
        m = predicted[predicted.treatment.==trt,:]
        d = dropmissing(observed[observed.treatment.==trt,:])
        sns.lineplot(m.tday, m.N, alpha=alpha, units=m.replicate, estimator=nothing, ax=ax[i])
        sns.lineplot(d.tday, d.count, ax=ax[i])
        if co_mean != nothing
            ax[i].axhline(co_mean)
        end
        ax[i].set(title=trt, xlabel="time (days)")
        if i==1
            ax[i].set_ylabel("# Daphnids")
        end
    end
    plt.tight_layout()
    sns.despine()
end

function plot_pdf(dist, range)
    p = [pdf(dist, x) for x in range]
    sns.lineplot(range, p)
    ax.set(xlabel="value", ylabel="probability density")
    ax.get_yaxis().set_visible(false)
    sns.despine()
end

"""
Plot population trajectories in different treatments, given predicted and observed data.
$(SIGNATURES)
"""
function nplot(predicted::DataFrame, observed::DataFrame; ci=95)
    
    trts = unique(predicted.treatment)
    fig, ax = plt.subplots(ncols=size(trts)[1], figsize=(20, 4), sharey=true)
    
    for (i, trt) in enumerate(trts)
        m = predicted[predicted.treatment.==trt,:]
        
        pest = by(m, :tday) do df
            est, lower, upper = get_estimate(df.N)
            DataFrame(
                est = est,
                lower = lower,
                upper = upper
            )
        end
        
        d = dropmissing(observed[observed.treatment.==trt,:])
        sns.scatterplot(d.tday, d.ind_per_L, ax=ax[i])
            
        sns.lineplot(pest.tday, pest.est, color="gray", estimator=nothing, ax=ax[i])
        sns.lineplot(pest.tday, pest.lower, color="lightgray", alpha=.5, ax=ax[i])
        sns.lineplot(pest.tday, pest.upper, color="lightgray", alpha=.5, ax=ax[i])
        
        ax[i].set(title=trt, xlabel="time (days)")
        
        if i==1
            ax[i].set_ylabel("# Daphnids")
        end
            
    end
        
    plt.tight_layout()
    sns.despine()
end

# functions for log(x+1)-transformed axes
logx = (py"lambda x:np.log(x+1)", py"lambda x:np.exp(x)-1")