using Random, StatsBase
using DataFrames
using DocStringExtensions
using CSV
include("Input.jl")
include("DEB.jl")
include("Inference.jl")

"""
Function to gather data from model object. Returns tidy DataFrame with a single row (= 1 time-step).
$(TYPEDSIGNATURES)
"""
function record_data(
        m::Model
    )
    
    # If there area any Fleas left, record individual-level data 
    if length(m.fleas)>0
        cum_repro_test_animal = [if f.unique_id==0 f.cum_repro end for f in m.fleas]
        carapace_length_test_animal = [if f.unique_id==0 f.L/f.shape_factor end for f in m.fleas]

        ind_statevar_mat = hcat([[
                f.food_taken .* m.timestep,
                f.S_A,
                f.S_C,
                f.cum_repro,
                f.L / f.shape_factor,
                f.functresp,
                f.Xi,
                f.D[1],
                f.D[2],
                f.D[3],
                f.f_stress,
                f.g_stress,
                f.k_M_stress_1,
                f.k_M_stress_2,
                f.k_J_stress,
                f.repro_stress_1,
                f.repro_stress_2
                ] for f in m.fleas]...)

        ingestion_rate_mean, 
        S_A_mean, 
        S_C_mean,
        cum_repro_mean,
        carapace_length_mean,
        functresp_mean,
        Xi_mean,
        D1_mean,
        D2_mean,
        D3_mean,
        f_stress_mean,
        g_stress_mean,
        k_M_stress_1_mean,
        k_M_stress_2_mean,
        k_J_stress_mean,
        repro_stress_1_mean,
        repro_stress_2_mean =
        [robust_mean(x) for x in eachrow(ind_statevar_mat)]

        # record life stage structure
        X_juv = []
        X_ad = []
        try
            push!(X_juv, sum([f.U_H < f.U_Hp for f in m.fleas]))
        catch
            push!(X_juv, 0)
        end
        try
            push!(X_ad, sum([f.U_H >= f.U_Hp for f in m.fleas]))
        catch
            push!(X_ad, 0)
        end    

        X_matrix = Array([])

        for j in 1:m.n_spec
            if j==1
                X_matrix = [get_species_abundance(m.fleas, j)]
            else
                X_matrix = hcat(X_matrix, get_species_abundance(m.fleas, j))
            end
        end 

        if ndims(X_matrix)==1
            X_matrix = add_dim(X_matrix)
        end

        model_out = hcat(
                    DataFrame(
                        tday            = [m.tday],
                        F               = [m.food_density],
                        cum_repro       = robust_mean(cum_repro_test_animal),
                        carapace_length = robust_mean(carapace_length_test_animal),
                        stressor_1_conc = [m.stressor_conc[1]],
                        stressor_2_conc = [m.stressor_conc[2]],
                        stressor_3_conc = [m.stressor_conc[3]],
                        ingestion_rate_mean = ingestion_rate_mean,
                        S_A_mean = S_A_mean,
                        S_C_mean = S_C_mean,
                        cum_repro_mean       = cum_repro_mean,
                        carapace_length_mean = carapace_length_mean,
                        functresp_mean = functresp_mean,
                        Xi_mean = Xi_mean,
                        D1_mean = D1_mean,
                        D2_mean = D2_mean,
                        D3_mean = D3_mean,
                        f_stress_mean = f_stress_mean,
                        g_stress_mean = g_stress_mean,
                        k_M_stress_1  = k_M_stress_1_mean,
                        k_M_stress_2  = k_M_stress_2_mean,
                        k_J_stress    = k_J_stress_mean,
                        repro_stress_1 = repro_stress_1_mean,
                        repro_stress_2 = repro_stress_2_mean,
                        X_juv = X_juv,
                        X_ad  = X_ad,
                        X     = [sum(X_matrix)]
                    ),
                    DataFrame(X_matrix)
                )

        return model_out
    else
        model_out = hcat(
              DataFrame(
                  tday            = [m.tday],
                  F               = [m.food_density],
                  cum_repro       = Array{Union{Float64,Missing}}([missing]),
                  carapace_length = Array{Union{Float64,Missing}}([missing]),
                  stressor_1_conc = [m.stressor_conc[1]],
                  stressor_2_conc = [m.stressor_conc[2]],
                  stressor_3_conc = [m.stressor_conc[3]],
                  ingestion_rate_mean = Array{Union{Float64,Missing}}([missing]),
                  S_A_mean = Array{Union{Float64,Missing}}([missing]),
                  S_C_mean = Array{Union{Float64,Missing}}([missing]),
                  cum_repro_mean       = Array{Union{Float64,Missing}}([missing]),
                  carapace_length_mean = Array{Union{Float64,Missing}}([missing]),
                  functresp_mean = Array{Union{Float64,Missing}}([missing]),
                  Xi_mean = Array{Union{Float64,Missing}}([missing]),
                  D1_mean = Array{Union{Float64,Missing}}([missing]),
                  D2_mean = Array{Union{Float64,Missing}}([missing]),
                  D3_mean = Array{Union{Float64,Missing}}([missing]),
                  f_stress_mean = Array{Union{Float64,Missing}}([missing]),
                  g_stress_mean = Array{Union{Float64,Missing}}([missing]),
                  k_M_stress_1  = Array{Union{Float64,Missing}}([missing]),
                  k_M_stress_2  = Array{Union{Float64,Missing}}([missing]),
                  k_J_stress    = Array{Union{Float64,Missing}}([missing]),
                  repro_stress_1 = Array{Union{Float64,Missing}}([missing]),
                  repro_stress_2 = Array{Union{Float64,Missing}}([missing]),
                  X_juv = [0.],
                  X_ad  = [0.],
                  X     = [0.]
              ),
              DataFrame(hcat(repeat([0.], length(m.fleas_dict_list))'...))
          )
    end
end

"""
Function to run the DEB-IBM, given global parameters in global_dict and individual-level parameters in fleas_dict.
$(SIGNATURES)
"""
function Run(
        global_dict::Dict{String,Any}, 
        fleas_dict::Dict{String,Any};
        cu_ts=nothing,
        ni_ts=nothing,
        zn_ts=nothing
    )
    
    #### ---- Initialize arrays to record data ---- ####
    tt = []
    Nt_juv = []
    Nt_ad = []
    Ft = []
    stressor_1_conc = []
    stressor_2_conc = []
    stressor_3_conc = []
    mean_damage_1 = []
    
    cum_repro = []
    carapace_length = []
    ingestion_rate = []
    functresp = []
    assim_rate = []
    
    max_steps = (global_dict["tmax"]*global_dict["timestep"])+1

    # if no exposure data has been passed on, use constant values given in global dictionary
    if cu_ts == nothing
        cu_ts = repeat([global_dict["stressor_conc"][1]], Integer(max_steps))
    end
    if ni_ts == nothing
        ni_ts = repeat([global_dict["stressor_conc"][2]], Integer(max_steps))
    end
    if zn_ts == nothing
        zn_ts = repeat([global_dict["stressor_conc"][3]], Integer(max_steps))
    end 
    
    # make a new model object
    m = Model(
            global_dict["remove_juveniles"],
            global_dict["renew_medium"],
            global_dict["frct_medium_renewal"],
            global_dict["timestep"],
            0, # tstep
            0, # tday
            global_dict["food_level"], # food_level
            global_dict["food_level"], # initial food_density
            global_dict["age_range_init"], 
            0, 
            [],
            fleas_dict,
            global_dict["stressor_conc"]
    )
    
    # generate initial juveniles
    for i in 1:global_dict["N0_juvenile"]
        NewFlea = make_flea(m, true, false)
        m.unique_id_count += 1
        push!(m.fleas, NewFlea)
    end
    # generate initial adults
    for i in 1:global_dict["N0_adult"]
        NewFlea = make_flea(m, false, true)
        m.unique_id_count +=1
        push!(m.fleas, NewFlea)
    end
    
    popsize=0
    
    # for every timestep
    for t in 1:max_steps
        # update time in days
        tday = (t-1)/m.timestep
        # read exposure concentrations
        m.stressor_conc[1] = cu_ts[t]
        m.stressor_conc[2] = ni_ts[t]
        m.stressor_conc[3] = zn_ts[t]
        # execute model step
        m = Step(m)
        # record population size
        popsize = length([f for f in m.fleas])
        
        if popsize > 5000
            @warn "Aborted simulation due to large population size"
            break
        end
        
        # collect time-series data
        try
            push!(Nt_juv, sum([f.U_H < f.U_Hp for f in m.fleas]))
        catch
            push!(Nt_juv, 0)
        end
        try
            push!(Nt_ad, sum([f.U_H >= f.U_Hp for f in m.fleas]))
        catch
            push!(Nt_ad, 0)
        end
        
        push!(Ft, m.food_density)
        push!(tt, tday)
        push!(stressor_1_conc, m.stressor_conc[1])
        push!(stressor_2_conc, m.stressor_conc[2])
        push!(stressor_3_conc, m.stressor_conc[3])
        try
            push!(cum_repro, mean([f.cum_repro for f in m.fleas]))
            push!(carapace_length, [if f.unique_id==0 f.L/f.shape_factor end for f in m.fleas][1])
            push!(ingestion_rate, mean([f.food_taken for f in m.fleas])*m.timestep)
            push!(functresp, mean(skipmissing([f.functresp for f in m.fleas])))
            push!(assim_rate, mean(skipmissing([f.S_A*m.timestep for f in m.fleas])))
        catch
            push!(cum_repro, missing)
            push!(carapace_length, missing)
            push!(ingestion_rate, missing)
            push!(functresp, missing)
            push!(assim_rate, missing)
        end
    end

    model_out = DataFrame(
            tday = tt, # time in days 
            N_juv = Nt_juv, # juveniles
            N_ad  = Nt_ad, # adults
            N = Nt_juv .+ Nt_ad, # population
            F = Ft, # food density
            stressor_1_conc = stressor_1_conc,
            stressor_2_conc = stressor_2_conc,
            stressor_3_conc = stressor_3_conc,
            cum_repro = cum_repro,
            carapace_length = carapace_length,
            ingestion_rate = ingestion_rate,
            functresp = functresp,
            assim_rate = assim_rate
        )
    
    model_out[isnothing.(model_out.carapace_length),:carapace_length] .= missing
    model_out[isnothing.(model_out.cum_repro),:cum_repro] .= missing
    
    # return time-series data as dataframe
    return model_out
end


"""
Function to simulate single-stressor toxicity tests. \n
$(SIGNATURES)
"""
function SingleStressorTest(
        global_dict::Dict{String, Any}, 
        fleas_dict::Dict{String, Any},
        stressor::Symbol,
        exposure;
        accepted_DEB = nothing,
        params_DEB = ["p_Am", "F_m", "y_EX", "mortality_constant"],
        accepted_stressor = nothing,
        replicates = 10
    )
    
    df = DataFrame()
    
    for rep in 1:replicates
        #### ---- Obtain parameter values ---- ####
        
        # make sure to always start with default parameters
        fleas_dict = copy(fleas_dict)
        
        # IF KDEs havee been passed on, sample parameter values from KDE
        # otherwise, use point estimates from fleas_dict
        
        # get DEB parameters
        if accepted_DEB != nothing
            sample = posterior_sample(accepted_DEB)
            for (i,par) in enumerate(params_DEB)
                    fleas_dict[par] = sample[i]
            end
            fleas_dict = fleas_dict |> convert_parameters |> estimate_egg_weight
        end
        
        # get toxicity parameters
        if accepted_stressor != nothing
            par_sample = posterior_sample(accepted_stressor)
            if stressor in [:Cu_ts, :Cu]
                fleas_dict = assign_params(fleas_dict, par_sample, string(stressor[1:2]))
            elseif stressor in [:Ni_ts, :Ni]
                fleas_dict = assign_params(fleas_dict, par_sample, string(stressor[1:2]))
            elseif stressor in [:Zn_ts, :Zn]
                fleas_dict = assign_params(fleas_dict, par_sample, string(stressor[1:2]))
            else
                error("Unknown stressor")
            end
        end

        # update compound parameters
        fleas_dict = fleas_dict |> convert_parameters |> estimate_egg_weight
        
        #### ---- Run the simulations ---- ####
        m = DataFrame()
        if typeof(exposure)==DataFrame
            for trt in unique(exposure.treatment)
                # get exposure time series
                stressor_ts = exposure[exposure.treatment.==trt,stressor]
                # run the model
                if stressor == :Cu_ts
                    m_trt = Run(global_dict, fleas_dict; cu_ts=stressor_ts)
                    m_trt.treatment = trt
                    m_trt.replicate = rep
                    append!(m, m_trt)
                elseif stressor == :Ni_ts
                    m_trt = Run(global_dict, fleas_dict; ni_ts=stressor_ts)
                    m_trt.treatment = trt
                    m_trt.replicate = rep
                    append!(m, m_trt)
                elseif stressor == :Zn_ts
                    m_trt = Run(global_dict, fleas_dict; zn_ts=stressor_ts)
                    m_trt.treatment = trt
                    m_trt.replicate = rep
                    append!(m, m_trt)
                else
                    error("Stressor has to be one of :Cu_ts, :Ni_ts, Zn_ts")
                end
            end
            
        elseif typeof(exposure)==Array{Float64,1}
            treatment_names = string(stressor) .* string.([x[1]-1 for x in collect(enumerate(exposure))])
            treatment_names[1] = "Co"
            # for every treatment
            for (i,trt) in enumerate(exposure)
                # run simulation with constant exposure
                if stressor == :Cu
                    global_dict["stressor_conc"][1] = trt
                    m_trt = Run(global_dict, fleas_dict)
                    m_trt[!,:treatment] .= treatment_names[i]
                    m_trt[!,:replicate] .= rep
                    append!(m, m_trt)
                elseif stressor == :Ni
                    global_dict["stressor_conc"][2] = trt
                    m_trt = Run(global_dict, fleas_dict)
                    m_trt[!,:treatment] .= treatment_names[i]
                    m_trt[!,:replicate] .= rep
                    append!(m, m_trt)
                elseif stressor == :Zn
                    global_dict["stressor_conc"][3] = trt
                    m_trt = Run(global_dict, fleas_dict)
                    m_trt[!,:treatment] .= treatment_names[i]
                    m_trt[!,:replicate] .= rep
                    append!(m, m_trt)
                else
                    error("Stressor has to be one of :Cu, :Ni, :Zn")
                end
            end
        else
            error("Exposure has to be of type DataFrame or Array{Float64, 1}, not "*string(typeof(exposure)))
        end
        # add a column for relative responses
        ref = by(m[m.treatment.=="Co",:], [:tday]) do df
            DataFrame(ref=mean(df.N))
        end
        
        m = join(m, ref, on=:tday)
        m.response = m.N ./ m.ref
        select!(m, Not(:ref))
    
        append!(df, m)

    end
    return df
end

"""
Function to simulate mixture treatments. 
$(TYPEDSIGNATURES)
"""
function MixtureTest(
        global_dict, 
        fleas_dict,
        exposure;
        accepted_DEB = nothing,
        params_DEB = ["p_Am", "F_m", "y_EX", "mortality_constant"],
        accepted_stressor_array = nothing,
        replicates = 10,
        treatments=["Co", "Mix1", "Mix2", "Mix3", "Mix4", "Mix5"]
    )
        
    df = DataFrame()
    for rep in 1:replicates
        m = DataFrame()
            
        # sample DEB parameters
        if accepted_DEB != nothing
            sample = posterior_sample(accepted_DEB)
            for (i,par) in enumerate(params_DEB)
                    fleas_dict[par] = sample[i]
            end
            fleas_dict = fleas_dict |> convert_parameters |> estimate_egg_weight 
        end
        
        # sample toxicity parameters
        if accepted_stressor_array!=nothing
            for (i,accepted_stressor) in enumerate(accepted_stressor_array)
                particle = posterior_sample(accepted_stressor)
                fleas_dict = assign_params(fleas_dict, particle, i)
            end
        end
        
        if typeof(exposure)==DataFrame
            for trt in unique(exposure.treatment)
                exp_trt = exposure[exposure.treatment.==trt,:]
                m_trt = Run(global_dict, fleas_dict; cu_ts=exp_trt[:,:Cu_ts], ni_ts=exp_trt[:,:Ni_ts], zn_ts=exp_trt[:,:Zn_ts])
                m_trt.treatment = trt
                append!(m, m_trt)
            end
        else
            for (i,trt) in enumerate(eachrow(exposure))
                global_dict["stressor_conc"][1] = trt[1]
                global_dict["stressor_conc"][2] = trt[2]
                global_dict["stressor_conc"][3] = trt[3]
                m_trt = Run(global_dict, fleas_dict)
                m_trt[!,:treatment] .= treatments[i]
                append!(m, m_trt)
            end
        end
        
        # Finally, add a column for relative responses
        ref = by(m[m.treatment.=="Co",:], [:tday]) do df
            DataFrame(ref=mean(df.N))
        end
            
        m = join(m, ref, on=:tday)
        m.response = m.N ./ m.ref
        select!(m, Not(:ref))
        m[!,:replicate] .= rep
        append!(df, m)
    end
    
    return df
end