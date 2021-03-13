using Random, Distributions
using DocStringExtensions

"
A mutable structure for models.
Contains model-level parameters and state variables. Intitialize with values from `global_dict`.
$(TYPEDFIELDS)
"
mutable struct Model
    remove_juveniles::Array{Int64,1} # weekdays on which juveniles are removed from the experiment
    renew_medium::Array{Int64,1} # weekdays on which medium is renewed
    frct_medium_renewal::Float64 # fraction of mefudium that is replaced at each renewal
    timestep::Int64
    tstep::Int64
    tday::Float64
    food_level::Float64
    food_density::Float64
    age_range_init::Int64
    unique_id_count::Int64
    fleas::Array
    fleas_dict::Dict
    stressor_conc::Array{Float64,1}
end

"""
A `toxicity_particle` contains all parameters related to toxic effects. \n
`k_d` is a one-dimenstional Array: There is only one dominant rate constant per stressor. \n
Hazard-related parameters are given in one-dimensional Arrays (one for each stressor). \n
Stress-related parameters are given in two dimensional Arrays (one for each stressor, one for each pMoA). \n
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct toxicity_particle
    k_d::Array{Float64,1}
    S_max::Array{Float64,2}
    S_ED50::Array{Float64,2}
    S_beta::Array{Float64,2}
    h_max::Array{Float64,1}
    h_ED50::Array{Float64,1}
    h_beta::Array{Float64,1}
end

"
A mutable structure for individuals. \n
Contains individual-level parameters and state variables. Initialize with values from `fleas_dict` using function `make_flea`.
$(TYPEDEF)
"
mutable struct Flea
    unique_id::Int64
    species::String
    cv::Float64
    
    # correction factors
    T_fact::Float64
    
    # scaled DEB parameters
    K::Float64
    J_XAm_rate::Float64
    U_Hb::Float64
    U_Hp::Float64
    k_M_rate::Float64
    k_J_rate::Float64
    g::Float64
    v_rate::Float64
    kap::Float64
    kap_R::Float64
    shape_factor::Float64
    mortality_constant::Float64
    starvation_hazard_rate::Float64
    sG::Float64
    h_a::Float64
    L0::Float64
    
    # embryonal parameters
    egg_weight::Float64
    L_embryo::Float64
    U_E_embryo::Float64
    
    # auxiliary
    time_between_molts::Float64

    crit_mass::Float64
    
    # toxicity parameters
    tox_particle::toxicity_particle
    
    # state variables
    die::Bool
    cause_of_death::String
    juvenile::Bool
    adult::Bool
    offspring_number::Int64
    cum_repro::Int64
    final_offspring::Int64
    molt_time::Float64
    develop_time::Float64
    U_E::Float64
    dU_E::Float64
    U_H::Float64
    dU_H::Float64
    U_R::Float64
    dU_R::Float64
    h_rate::Float64
    dh_rate::Float64
    q_acceleration::Float64
    dq_acceleration::Float64
    L::Float64
    dL::Float64
    L_max::Float64
    Xi::Float64
    functresp::Float64
    food_taken::Float64
    S_A::Float64
    S_C::Float64
    e_scaled::Float64
    
    # _int values for everything that can be changed by a stressor
    J_XAm_rate_int::Float64
    k_J_rate_int::Float64
    g_int::Float64
    k_M_rate_int::Float64
    kap_R_int::Float64
    
    # toxicity-related state variables per stressor
    D::Array{Float64,1} # damage
    h_z::Array{Float64,1} # hazard rate
    p_z::Array{Float64,1} # survival probability 
    stress::Array{Float64,2} # stress
    
    # toxicity-related state variables per pMoA
    f_stress::Float64
    g_stress::Float64
    k_M_stress_1::Float64
    k_M_stress_2::Float64
    k_J_stress::Float64
    repro_stress_1::Float64
    repro_stress_2::Float64
    
    # individual modifier
    scatter_multiplier::Float64   
end

function make_tox_particle(tox_particle_dict::Dict)
    tox_particle = toxicity_particle(
    tox_particle_dict["k_d"],
    tox_particle_dict["S_max"],
    tox_particle_dict["S_ED50"],
    tox_particle_dict["S_beta"],
    tox_particle_dict["h_max"],
    tox_particle_dict["h_ED50"],
    tox_particle_dict["h_beta"]
    )
    return tox_particle
end

"
A function to initialize a Flea. \n
Choose start `juvenile=1` or start `adult=1` for initial individuals.
$(SIGNATURES)
"
function make_flea(
        m::Model,
        start_juvenile::Bool,
        start_adult::Bool
    )
    
    fleas_dict = m.fleas_dict
    
    # generate Flea object
    f = Flea(
        m.unique_id_count,
        fleas_dict["species"], 
        fleas_dict["cv"],
        # correction factors
        1.,
        
        # individual-level parameters
        fleas_dict["K"],
        fleas_dict["J_XAm_rate"],
        fleas_dict["U_Hb"],
        fleas_dict["U_Hp"],
        fleas_dict["k_M_rate"],
        fleas_dict["k_J_rate"],
        fleas_dict["g"],
        fleas_dict["v_rate"],
        fleas_dict["kap"],
        fleas_dict["kap_R"],
        fleas_dict["shape_factor"],
        fleas_dict["mortality_constant"],
        fleas_dict["starvation_hazard_rate"],
        fleas_dict["sG"],
        fleas_dict["h_a"],
        fleas_dict["L0"],
        
        # embryonal parameters
        fleas_dict["egg_weight"],
        fleas_dict["L_embryo"],
        fleas_dict["U_E_embryo"],
        
        # other
        fleas_dict["time_between_molts"],
        fleas_dict["crit_mass"],
        
        # toxicity parameters
        make_tox_particle(fleas_dict["toxicity_particle"]),
        
        # state variables
        false, # die?
        "None",
        0, # juvenile 
        0, # adults
        0, # offspring_number
        0, # final_offspring
        0, # cum_repro
        0., # molt_time
        0., # develop_time
        0., # U_E
        0., # dU_E
        0., # U_H
        0., # dU_H
        0., # U_R
        0., # dU_R
        0., # h_rate
        0., # dh_rate
        0., # q_acceleration
        0., # dq_acceleration
        0., # L
        0., # dL
        0., # L_max
        0., # Xi
        0., # functresp
        0., # food_taken
        0., # S_A
        0., # S_C
        0., # e_scaled
        
        # _int values
        fleas_dict["J_XAm_rate"],
        fleas_dict["k_J_rate"],
        fleas_dict["g"],
        fleas_dict["k_M_rate"],
        fleas_dict["kap_R"],
        
        # toxicity-related state variables per stressor
        [0., 0., 0.], # damage
        [0., 0., 0.], # hazard rate
        [1., 1., 1.], # survival probability
        repeat([0., 0., 1., 0.], 1, 3), # stress
        
        # toxicity-related state variables per pMoA
        1., # f_stress
        0., # g_stress
        0., # k_M_stress_1
        0., # k_M_stress_2
        0., # k_J_stress
        0., # repro_stress_1
        0., # repro_stress_2
        
        # individual variability
        0. # scatter_multiplier
    )
    # TODO at some point: implement temperature correction
    f.T_fact = 1
    
    f = individual_variability(f)
    
    if start_juvenile
        f = init_random_age(f, m)
    elseif start_adult
        f = init_adulthood(f, m)
    end
    
    return f
end

""" 
Initialize a juvenile with age <= age _ range _ init. Age is given in units of model timesteps.
Initial juveniles are typically 0-24 hours old.
$(SIGNATURES)
"""
function init_random_age(f::Flea, m::Model)
    f.develop_time = -1
    f.juvenile = false
    
    age = rand(1:m.age_range_init)
    
    for a in 1:age
        f.functresp = 1
        f, m = feed(f, m, false)
        f = dU_E(f)
        f = dU_H(f)
        f = dU_R(f)
        f = dL(f)
        f = update(f, f.cv, m.timestep)
    end
    return f
end

"
Simulate the live of an adult before the start of an experiment.
It is assumed that the individual could always feed ad libitum, therefore the scaled functional response
is fixed to 1. 
Investment in reproduction is accounted for and individuals have the potential to produce offspring immediately after the onset of the experiment.

Usage:
```
f = init_adulthood(f, m)
```
"
function init_adulthood(f::Flea, m::Model)
    init_time = 0
    while (f.shape_factor*10<3) & (init_time < m.timestep*20)
        init_time +=1
        f,m = feed(f, m, false)
        f = dU_E(f) 
        f = dU_H(f)
        f = dU_R(f)
        f = dL(f)
        f = update(f, f.cv, m.timestep)
        
        if f.adult
            f.molt_time = f.molt_time + (f.T_fact/m.timestep)
            
            # emulate molting + reproduction
            if f.molt_time > f.time_between_molts
                f.molt_time = 0
                f.final_offspring = Int64(floor((f.U_R/f.kap_R)/(f.egg_weight/f.T_fact)))
                
                # TODO:
                # this line is not included in netlogo model
                # but it makes sense?
                # think about whether to include this
                f.U_R = f.U_R - ((floor(f.U_R/(f.egg_weight/f.T_fact))) * (f.egg_weight/f.T_fact))
            end
        end
        
    # TODO
    # figure out whether this is really necessary (relict of netlogo model)
    # annotate and remove magic numbers
    f.molt_time =  -1*rand(1:4)
    f.U_R = (f.egg_weight/f.T_fact) * (15 * f.molt_time) / f.time_between_molts
    
    end
    return f
end

"""
Feeding of individuals.
Returns updated model-level (food density) and individual-level (conditional food density, functional response) state variables.
$(SIGNATURES)
"""
function feed(;
        U_H::Float64, 
        U_Hb::Float64, 
        K::Float64, 
        L::Float64, 
        T_fact::Float64, 
        J_XAm_rate::Float64,
        f_stress::Float64,
        timestep::Int64, 
        food_density::Float64, 
        update_food_density::Bool
    )
    # Embryos don't feed
    if U_H <= U_Hb
        functresp = 0
        Xi = 0
    # Juveniles and adults do
    else
        Xi = food_density
        # functional response
        functresp = Xi/(K + Xi)
    end
    
    # to convert scaled functional response to amount of ingested food
    # - multiply with surface area: the result is identical with assimilation flux S_A,
    # - multiply with the surface area-specific maximum ingestion rate
    # - multiply with temperature correction factor
    food_taken = (functresp * L^2 * T_fact * J_XAm_rate) / timestep

    if update_food_density
        # subtract ingested food from environmental food density
        # max() is applied just to make sure that floating point errors don't mess things up 
        # (if food_density ever becomes negative, the error will propagate almost irreversibly)
        food_density = max(0, food_density - food_taken)
    end
    
    return Xi, functresp, food_taken, food_density 
end

function feed(f::Flea, m::Model, update_food_density::Bool)
    f.Xi, f.functresp, f.food_taken, m.food_density = feed(
            U_H          = f.U_H, 
            U_Hb         = f.U_Hb, 
            K            = f.K, 
            L            = f.L, 
            T_fact       = f.T_fact,
            J_XAm_rate   = f.J_XAm_rate,
            f_stress     = f.f_stress,
            timestep     = m.timestep, 
            food_density = m.food_density, 
            update_food_density = update_food_density)
    return f, m
end


"""
Induction of individual variability, based on scatter_multiplier.
"""
function individual_variability(f::Flea)
    f.scatter_multiplier = exp(rand(Normal(0, f.cv)))
    
    f.g    = f.g / f.scatter_multiplier
    f.U_Hb = f.U_Hb / f.scatter_multiplier
    f.U_Hp = f.U_Hp / f.scatter_multiplier
    f.K    = f.K * f.scatter_multiplier
    f.J_XAm_rate = f.J_XAm_rate * f.scatter_multiplier
    f.L   = f.L_embryo * f.scatter_multiplier
    f.U_E = f.U_E_embryo * f.scatter_multiplier
    f.juvenile = true
    f.adult    = false
    f.offspring_number = 0
    f.U_R  = 0
    f.dU_R = 0
    f.U_H  = f.U_Hb
    f.dU_H = 0
    f.h_rate  = 0
    f.dh_rate = 0
    f.q_acceleration  = 0
    f.dq_acceleration = 0
    
    # store _int values for everything that can be changed by a stressor
    f.k_J_rate_int = f.k_J_rate
    f.g_int        = f.g
    f.k_M_rate_int = f.k_M_rate
    f.kap_R_int    = f.kap_R
    
    return f
end

"""
Calculate change in scaled maturity.
"""
function dU_H(f::Flea)    
    if f.U_H < f.U_Hp/f.T_fact
        f.dU_H = (1-f.kap) * f.S_C - (f.T_fact * f.k_J_rate) * f.U_H
    else    
        f.dU_H = 0
    end
    return f
end

"""
Calculate change in scaled reserves.
Also: Assimilation flux, mobilization flux, scaled reserve density.
"""
function dU_E(f::Flea)
    if f.juvenile
        f.S_A      = f.functresp * f.L^2 * f.f_stress
        f.e_scaled = f.T_fact * f.v_rate * (f.U_E/(f.L^3))
        f.S_C      = f.L^2 * (f.g*f.e_scaled/(f.g+f.e_scaled)) * (1+(f.L/(f.g*(f.T_fact*f.v_rate/(f.g*f.T_fact*f.k_M_rate)))))
        f.dU_E     = f.S_A - f.S_C
    end
    return f
end

"""
Calculate change in reproduction buffer
"""
function dU_R(f::Flea)
    if f.U_H > f.U_Hp
        f.dU_R = (1-f.kap) * f.S_C - (f.T_fact * f.k_J_rate) * (f.U_Hp / f.T_fact)
    else
        f.dU_R = 0
    end
    return f
end

"""
Calculate change in scaled size.
"""
function dL(f::Flea)
    
    f.dL = (1/3) * 
    (((f.T_fact * f.v_rate)/(f.g*(f.L^2)) * f.S_C) - (
            f.T_fact * f.k_M_rate) * f.L)

    # starvation if growth is negative
    """
    if f.e_scaled < f.L / (f.v_rate / (f.g * f.k_M_rate))
        f.dL = 0
        
        # juveniles divert maturation energy into maintenance
        if f.dU_H < (f.U_Hp / f.T_fact)
            f.dU_H = (1-f.kap) * f.e_scaled * f.L^2 - f.k_J_rate * f.U_Hp - f.kap * f.L^2 * (f.L / (f.v_rate / (f.g * f.k_M_rate)) - f.e_scaled)
        # adults divert reproduction energy into maintenance
        else
            f.dU_R = (1-f.kap) * f.e_scaled * f.L^2 - f.k_J_rate * f.U_Hp - f.kap * f.L^2 * (f.L / (v.v_rate / (f.g * f.k_M_rate)) - f.e_scaled)
        end
        
        f.dU_E = f.S_A - f.e_scaled * f.L ^ 2
        
        if f.U_H < (f.U_Hp / f.T_fact)
            if f.dU_H < 0
                f.die = true
                end
        else
            if f.dU_R < 0
                f.die = true
            end
        end
    end
    """
    
    return f
end

"""
Calculate ageing acceleration
"""
function dq_acceleration(f::Flea)
    f.dq_acceleration = (f.q_acceleration * ((f.L^3)/(((f.T_fact*f.v_rate)/(f.g*(f.T_fact*f.k_M_rate)))^3))*f.sG+f.h_a)*f.e_scaled*((f.T_fact*f.v_rate/f.L)-((3/f.L)*f.dL))-((3/f.L)*f.dL)*f.q_acceleration
    return f
end

function dh_rate(f::Flea)
    f.dh_rate = f.q_acceleration - ((3/f.L)*f.dL)*f.h_rate
    return f
end

function update(f::Flea, cv::Float64, timestep::Int64)
    # update main state variables
    f.U_E = f.U_E + f.dU_E/timestep
    f.U_H = f.U_H + f.dU_H/timestep
    f.U_R = max(f.U_R + f.dU_R/timestep, 0)
      f.L = max(0, f.L + f.dL/timestep)
    f.L_max = max(f.L, f.L_max)
    # update ageing and hazard rate
    if f.U_H > f.U_Hb
        f.q_acceleration = f.q_acceleration + f.dq_acceleration/timestep
        f.h_rate = f.h_rate + f.dh_rate/timestep
    end
    # transition from embryo to juvenile
    if f.U_H >= f.U_Hb/f.T_fact
        f.juvenile = true
    end
    # transition from juvenile to adult
    if f.U_H > f.U_Hp/f.T_fact
        f.adult = true
    end
    # check development status
    if f.juvenile==false
        f.develop_time = f.develop_time - f.T_fact
        # introduce individual variability in embryos
        if f.develop_time <= 0
            f = individual_variability(f)
        end
    end
    return f
end

function reset(f::Flea)
    f.adult    = 0
    f.juvenile = 0
    f.offspring_number = 0
    f.L    = 0
    f.U_E  = 0
    f.U_H  = 0
    f.U_R  = 0
    f.dU_R = 0
    f.h_rate  = 0
    f.dh_rate = 0
    f.q_acceleration  = 0
    f.dq_acceleration = 0
    f.S_A = 0
    f.S_C = 0
    f.dL  = 0
    f.L_max = 0 # maximum structural size
    f.develop_time = f.T_fact - 1
    f.die = false
    return f
end

"""
Inheritance attributes.
$(SIGNATURES)
"""
function inherit(parent::Flea, child::Flea)
    child.g    = parent.g
    child.U_Hb = parent.U_Hb
    child.U_Hp = parent.U_Hp
    child.K    = parent.K
    child.J_XAm_rate = parent.J_XAm_rate
    child.L_embryo   = parent.L_embryo
    child.U_E_embryo = parent.U_E_embryo
    
    return child
end

function reproduce(f::Flea, m::Model)
    # molting controls timing of reproduction
    f.molt_time = f.molt_time + (f.T_fact/m.timestep)
    if f.molt_time > f.time_between_molts
        # determine no. of eggs based on reproduction buffer and reproduction efficiency
        num_eggs = (f.U_R*f.kap_R)/(f.egg_weight/f.T_fact)
        # convert to integer
        if isnan(num_eggs)
            f.offspring_number = 0
        else
            f.offspring_number = Int64(floor(num_eggs))
            f.cum_repro = f.cum_repro + Int64(floor(num_eggs))
        end
        # reset molting
        f.molt_time = 0
        # the subtracted amount of reserve does not take kap_R into account 
        # (energy is spent, but not effectively turned into offspring)
        f.U_R = f.U_R - ((floor(f.U_R/(f.egg_weight/f.T_fact))) * (f.egg_weight/f.T_fact))
        # create new offspring
        for i in 1:f.offspring_number
            NewFlea = make_flea(m, false, false)
            NewFlea = inherit(f, NewFlea)
            NewFlea = reset(NewFlea)
            NewFlea = individual_variability(NewFlea)
            m.unique_id_count += 1
            # add new flea to fleas array
            push!(m.fleas, NewFlea)
        end
    end
    return f, m
end

function death(f::Flea, m::Model)
    # death by starvation
    condition_1a = (f.L^3) < f.crit_mass * (f.L_max^3)
    condition_1b = rand() < 1 - ((1-f.starvation_hazard_rate)^(f.T_fact/m.timestep))
    
    if condition_1a & condition_1b
        f.cause_of_death = "starvation"
        f.die = true
    end
    
    # death by ageing
    f.h_rate = min(f.h_rate, 1)
    condition_2 = rand() < 1-((1-f.h_rate)^(1/m.timestep))
        
    if condition_2
        f.cause_of_death = "age"
        f.die = true
    end
    
    # experimental removal of animals
    condition_3a = true in [floor(m.tday)%7==x for x in m.remove_juveniles]
    condition_3b = f.unique_id != 0
    
    if condition_3a & condition_3b
        f.cause_of_death = "removal"
        f.die = true
    end
    
    return f
end

"""
Hazard rate dependent on scaled reserve density.
$(SIGNATURES)
"""
function density_dependent_mortality(f::Flea, m::Model)
    if (f.juvenile==true) & (f.adult==false)
        mortality_probability = (1 - f.e_scaled) * (1-(1-f.mortality_constant)^(1/m.timestep))
        if rand() < mortality_probability
            f.cause_of_death = "density"
            f.die = true
        end
    end
    return f
end

"""
change in damage
$(SIGNATURES)
"""
function dD(D, k_d, stressor_conc)
    return D + k_d * (stressor_conc - D)
end

"""
log-logistic response
$(SIGNATURES)
"""
function log_logistic(x, upper_limit, EC50, slope)
    return upper_limit * (1/(1+(x/EC50)^slope))
end

"""
Lethal and sublethal effects of an arbitrary number of stressors. All stressors can have all pMoAs simultaneously.
$(SIGNATURES)
"""
function stressor_toxicity(f::Flea, m::Model)
    
    #### Step 1: Reset physiological stress ####
    f.f_stress = 1
    f.g_stress = 0
    f.k_M_stress_1 = 0
    f.k_J_stress = 0
    f.repro_stress_1 = 0
    f.repro_stress_2 = 0
    
    #### Step 2: update scaled damage ####
    f.D = max.(0, dD.(f.D, f.tox_particle.k_d, m.stressor_conc))
    
    #### Step 3: lethal effects ####
    ## Step 3.1: calculate hazard rate ##
    f.h_z = log_logistic.(
        f.D, 
        f.tox_particle.h_max, 
        f.tox_particle.h_ED50, 
        .-f.tox_particle.h_beta
    )
    
    ## Step 3.2: convert to survival probability ##
    f.p_z = exp.(.- f.h_z ./m.timestep)
    
    ## Step 3.3: stochastic death ##
    dice = rand(length(f.p_z))
    if sum(dice .> f.p_z) .>=1
        f.die=true
        f.cause_of_death="stressor"
    end
    
    # apply log-logistic function for every stressor
    f.stress = log_logistic.(
        f.D', 
        f.tox_particle.S_max, 
        f.tox_particle.S_ED50, 
        .-f.tox_particle.S_beta
    )
    
    #### Step 4: calculate mixture stress for every pMoA ####
    # Growth
    f.g_stress     = prod(1 .+ f.stress[1,:])
    f.k_M_stress_2 = prod(1 ./ (1 .+ f.stress[1,:]))

    # Maintenance
    f.k_J_stress   = prod(1 .+ f.stress[2,:])
    f.k_M_stress_1 = prod(1 .+ f.stress[2,:])

    # Assimilation
    f.f_stress = prod(1 .- f.stress[3,:])

    # Reproduction
    f.repro_stress_1 = prod(1 ./ (1 .+ f.stress[4,:]))
    
    #### Step 5: apply stress to physiological parameters ####
    # note that f_stress ins applied in `dU_E()`
    f.g = f.g_int * f.g_stress
    f.k_M_rate = f.k_M_rate_int * f.k_M_stress_1 * f.k_M_stress_2
    f.k_J_rate = f.k_J_rate_int * f.k_J_stress
    f.kap_R = f.kap_R_int * f.repro_stress_1
    
    return f
end

"""
Update environmental state variables. \n
Includes renewal of medium (=removal of food) and addition of food.
$(SIGNATURES)
"""
function update_environment(m::Model)
    # add food daily
    if m.tstep % m.timestep == 0
        # renew media
        if true in [floor(m.tday)%7==x for x in m.renew_medium]
             m.food_density = m.food_density * (1 - m.frct_medium_renewal)
        end            
        m.food_density = m.food_density + m.food_level * exp(rand(Normal(0.0, 0.2)))
    end
    return m
end


"""
Execute a single model step.
"""
function Step(m::Model)
    # randomize order in which individuals are called
    shuffle!(m.fleas)
    # execute individual schedule
    for f in m.fleas
        
        f = stressor_toxicity(f, m)
        f, m = feed(f, m, true)
        if f.juvenile
            f = f |> 
            dU_E |>
            dU_H |> 
            dU_R |> 
            dL |> 
            dq_acceleration |> 
            dh_rate
        end
        f = density_dependent_mortality(f, m)
        f = update(f, f.cv, m.timestep)
        if f.adult
            f, m = reproduce(f, m)
        end
        f = death(f, m)
    end
    # add food, renew media etc
    m = update_environment(m)
    # remove dead individuals
    filter!(f->f.die==false, m.fleas)
    m.tstep +=1
    m.tday = m.tstep/m.timestep
    return m
end
