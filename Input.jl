#Input.jl
"""
Handling input for DEB-IBM simulations.
Contains functions for parameter conversions and corrections, as well as dictionaries with AMP values and default simulation settings.
"""

using DocStringExtensions

#### ----- default global dictionaries ---- ####

# parameters for population experiment with D. magna conducted in early 2019
population_dict = Dict(
    "remove_juveniles" => [], # weekdays on which to remove juveniles
    "renew_medium" => [0, 3], # weekdays on which to remove medium
    "frct_medium_renewal" => 0.25, # fraction of medium to replace per renewal
    "N0_juvenile" => 8, # initial juveniles
    "N0_adult"    => 2, # initial adults
    "timestep"    => 24, # number of timesteps per day
    "tmax"        => 42, # max. simulated time (days)
    "food_level"  => 8e7, # food level (# cells/L/day)
    "age_range_init" => 24, # maximum initial age of juveniles (timesteps)
    "stressor_conc" => [0., 0., 0.]
)

# settings of a life-table experiment as conducted in late 2019/early 2020
life_table_dict = Dict(
    "remove_juveniles" => [0, 2, 5], # remove juveniles 3x per week
    "renew_medium"     => [0, 2, 5], # renew media on the same days
    "frct_medium_renewal" => 1.00, # 100% of media are renewed
    "N0_juvenile" => 1, # simulation starts with 1 juvenila
    "N0_adult"    => 0, 
    "timestep"    => 24, # number of timesteps per day
    "tmax"        => 21, # test lasts for 21 days
    "food_level"  => 5., # food level in mgC//L3/day
    "age_range_init" => 24, # maximum initial age of juveniles (timesteps)
    "stressor_conc" => [0., 0., 0.]
)

# reference scenario for ecorelevance project
ecorelevance_reference_dict = Dict(
    "remove_juveniles" => [],
    "renew_medium" => [0,3], # renew media on mondays and thursdays
    "frct_medium_renewal" => 0.25, # renew 25% media
    "N0_juvenile" => 3, # start with 3 juveniles
    "N0_adult" => 1, # and one adult
    "age_range_init" => 48,
    "timestep" => 24,
    "tmax" => 56,
    "food_level" => 2.5/1000, # food level in mgC/cm^3/day
    "stressor_1_conc" => 0, # concentration of stressor 1
    "stressor_2_conc" => 0, # concentration of stressor 2
    "stressor_3_conc" => 0  # concentration of stressor 3
)


"""
Conversion of AMP parameters for use in IBM.
See user manual of the generic DEB-IBM implementation for explanations.
$(SIGNATURES)
"""
function convert_parameters(fld)
    # maintenance rate follows from volume-specific somatic structure
    # and volume-specific costs for structure
    fld["k_M_rate"]    = fld["p_M"]  / fld["E_G"]
     # compound parameters following from scaling by assimilation rate
    fld["g"]           = fld["E_G"]  * fld["v_rate"] / (fld["kap"] * fld["p_Am"])
    fld["U_Hb"]        = fld["E_Hb"] / fld["p_Am"]
    fld["U_Hp"]        = fld["E_Hp"] / fld["p_Am"]
    # functional response parameters derived from primary parameters
    fld["J_XAm_rate"]  = fld["p_Am"] / fld["y_EX"]
    fld["K"]           = fld["J_XAm_rate"] / fld["F_m"]
    return fld
end

"""
Estimate embryonal parameters using bisection method.
"""
function estimate_egg_weight(
        timestep::Int64,
        v_rate::Float64,
        U_Hb::Float64,
        kap::Float64,
        g::Float64,
        k_J_rate::Float64,
        k_M_rate::Float64
)
    embryo_timestep = timestep * 1000
    lower_bound = 0
    upper_bound = 1
    maxiter = 100
    egg_weight = 0.5
    U_H_embryo = 0
    U_E_embryo = 0
    e_scaled_embryo = 0
    e_ref = 1
    L_embryo = 0.00001
    
    for i in 1:maxiter
        egg_weight = .5 * (lower_bound + upper_bound)
        L_embryo = 0.00001
        U_E_embryo = egg_weight
        U_H_embryo = 0
        e_scaled_embryo = v_rate * (U_E_embryo/L_embryo)^3
        
        while (U_H_embryo < U_Hb) & (e_scaled_embryo >1)
             # calculate reserve density based on energy conductance, scaled reserves and structural volume
            e_scaled_embryo = v_rate * (U_E_embryo/L_embryo^3)
            # calculate corresponding mobilization flux
            S_C_embryo      = L_embryo^2 * (g * e_scaled_embryo / ( g + e_scaled_embryo)) * (1 + (L_embryo / (g * (v_rate / (g * k_M_rate)))))
            # calculate change in state varaibles
            dU_E_embryo     = -1 * S_C_embryo
            dU_H_embryo     = (1-kap)*S_C_embryo - k_J_rate*U_H_embryo
            dL_embryo       = ((1/3) * (((v_rate / (g*(L_embryo^2))) * S_C_embryo) - k_M_rate*L_embryo))
            # update state variables
            U_E_embryo = U_E_embryo + dU_E_embryo / embryo_timestep
            U_H_embryo = U_H_embryo + dU_H_embryo / embryo_timestep
            L_embryo   = L_embryo   + dL_embryo   / embryo_timestep
            e_scaled_embryo = v_rate * (U_E_embryo / L_embryo^3)        
        end
        
        # in order to accept an estimate:
        # 1. estimated reserve density of embryo has to be close to one 
        # (to fulfil the criterion that embryos are mostly reserve)
        # 2. estimated maturity has to be larger or equal to maturity at birth (by definition)
        if (e_scaled_embryo < 0.01+e_ref) & (e_scaled_embryo > -0.01+e_ref) & (U_H_embryo >= U_Hb)
            @goto escape_label # jump to @label escape_label
        end
        # if estimate has not been accpeted:
        # adjust upper/lower bound based on result of previous estimation
        if U_H_embryo > U_Hb
            upper_bound = egg_weight
        else
            lower_bound = egg_weight
        end
        #if i == maxiter
        #    error("Embryo submodel is not converging. Smaller timesteps could help.")
        #end
    end
    @label escape_label
    
    return egg_weight, L_embryo, U_E_embryo, U_H_embryo
end

"""
Estimation of intitial state variables of an Embryo, determining the amount of energy
that has to be invested per offspring.
$(SIGNATURES)
"""
function estimate_egg_weight(fleas_dict::Dict)
    # estimate embryonal parameters
    fleas_dict["egg_weight"], 
    fleas_dict["L_embryo"], 
    fleas_dict["U_E_embryo"],
    fleas_dict["U_H_embryo"] = estimate_egg_weight(
        24, 
        fleas_dict["v_rate"], 
        fleas_dict["U_Hb"],
        fleas_dict["kap"],
        fleas_dict["g"],
        fleas_dict["k_J_rate"],
        fleas_dict["k_M_rate"])
    return fleas_dict
end

#### ---- default flea dictionaries ---- ####
"""
Parameter values for Daphnia magna according to AMP database.
Last checked: Jan 2020
Toxicity parameters from Daphnia magna models (Pereira et al. 2019; Karel Vlaeminck, unpublished).
"""
DM_AMP = Dict(
    "species" => "DM_AMP",
    "cv"      => 0.1,
    
    # AMP parameters
    "p_Am" => 313.169, #
    "F_m"  => 30., # filtration rate (cm^3/day*cm^2)
    "E_Hb"       => 0.05464, 
    "E_Hp"       => 1.09,
    "p_M"        => 1200,
    "E_G"        => 4400,
    "k_J_rate"   => 0.2537,
    "v_rate"     => 0.1858,
    "kap"        => 0.5809,
    "kap_R"      => 0.95,
    "shape_factor"       => 0.264,
    
    "sG"  => -0.3,
    "h_a" => 0.0002794, 
    "time_between_molts" => 2.324,
   
    # embryonal parameters
    "L0"  => 0.0001,
    "egg_weight" => 3.771781921386719E-4,
    "L_embryo"   => 0.023340451224020977,
    "U_E_embryo" => 6.039673323323045E-5,
    
    "starvation_hazard_rate" => 0.,
    "crit_mass" => 0.,
  
    # Toxicity parameters
    # Toxicity parameters
    "pmoas"  => ["Growth", "Maintenance", "Assimilation", "Reproduction"],
    "k_d"    => [ 0., 0., 0.],
    "h_max"  => [ 0., 0., 0.],
    "h_beta" => [ 0., 0., 0.],
    "h_ED50" => [ 0., 0., 0.],
    "S_max"  => [ 0., 0., 0.],
    "S_beta" => [ 0., 0., 0.],
    "S_ED50" => [ 0., 0., 0.] 
)

"""
Parameter values for Daphnia magna, fed with mixture of Pseudokirchneriella subcapitata and Chlamydomonas reinhardtii, simulated as a single food source.
Food densities in # cells/L.
"""
DM1_2019 = Dict(
    "species" => "DM1_2019",
    "cv"      => 0.1,
    
    # corrected parameters
    
    "p_Am" => 260.4, # (CI = (260.3, 309.5))
    "F_m"  => 10.0, # (CI = (8, 11))
    "y_EX" => 3.7e-6, # (CI = (3e-6, 4e-6))
    "mortality_constant" => 0.41, # (CI = (0.3, 0.63))
    
    # AMP parameters
    
    "E_Hb"       => 0.05464, 
    "E_Hp"       => 1.09,
    "p_M"        => 1200,
    "E_G"        => 4400,
    "k_J_rate"   => 0.2537,
    "v_rate"     => 0.1858,
    "kap"        => 0.5809,
    "kap_R"      => 0.95,
    "shape_factor"       => 0.264,
    
    "sG"  => -0.3,
    "h_a" => 0.0002794, 
    "time_between_molts" => 2.324,
   
    # embryonal parameters
    "L0"  => 0.0001,
    "egg_weight" => 3.771781921386719E-4,
    "L_embryo"   => 0.023340451224020977,
    "U_E_embryo" => 6.039673323323045E-5,
    
    "starvation_hazard_rate" => 0.,
    "crit_mass" => 0.,
  
    "toxicity_particle" => Dict(
        "stressors" => ["Cu", "Ni", "Zn"],
        "pMoAs" => ["Growth", "Maintenance", "Assimilation", "Reproduction"],
        # one k_d for Cu, Ni and Zn
        "k_d" => [0., 0., 0.],
        # one parameter for every stressor (rows) and every pMoA (columns)
        "S_max" => hcat(
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.]
        ),
        "S_ED50" => fill(Inf, 4, 3),
        "S_beta" => fill(2., 4, 3),
        "h_max" => [0., 0., 0.],
        "h_ED50" => [Inf, Inf, Inf],
        "h_beta" => [2., 2., 2.]
    )
)


DL_AMP = Dict(
    "species" => "DL_AMP",
    "cv"      => 0.1,
    
    # AMP parameters
    "E_Hb"       => 0.01746, 
    "E_Hp"       => 1.401,
    "p_M"        => 5686.83,
    "E_G"        => 4400.59,
    "k_J_rate"   => 0.002,
    "v_rate"     => 0.020455,
    "kap"        => 0.5568,
    "kap_R"      => 0.95,
    "shape_factor"       => 0.18946,
    "mortality_constant" => 0,
    "crit_mass"          => 0.4,
    "starvation_hazard_rate" => 0.9,
    "sG"  => 0.001,
    "h_a" => 8.603e-5,
    "L0"  => 0.0001,  

    # food-dependent parameters
    "p_Am" => 558.735, # J/d*cm^2
    "y_EX" => 1e-3, # yield of reserve (J) on food (# cells/cm^3)
    "F_m"  => 270.833, # filtration rate (cm^3/day*cm^2)
    
    # embryonal parameters
    "egg_weight" => 3.771781921386719E-4,
    "L_embryo"   => 0.023340451224020977,
    "U_E_embryo" => 6.039673323323045E-5,
    
    # reproductive timing
    "time_between_molts" => 2.324,
    
    # Toxicity parameters
    # compound 1
    "compound_1" => "Copper", 
    "pmoa_1"     => "Growth", 
    "h_b_1"      => 0/24, # background hazard rate
    "k_d_1"      => 0, # dominant rate constant
    "k_k_1"      => 0, # ?
    "m_n_1"      => 0, # ?
    "s_s_1"      => 0,
    "s_d_1"      => 0, # ?
    "A_eq_1" => 0,
    "B_eq_1" => 0,
    "C_eq_1" => 0,
    
    # compound 2
    "compound_2" => "Nickel",
    "pmoa_2"     => "Maintenance",
    "h_b_2"      => 0,
    "k_d_2"      => 0,
    "k_k_2"      => 0,
    "m_n_2"      => 0,
    "s_s_2"      => 0,
    "s_d_2"      => 0,
    "A_eq_2"     => 0,
    "B_eq_2"     => 0,
    "C_eq_2"     => 0,
    
    # compound 3
    "compound_3" => "Zinc",
    "pmoa_3"     => "Assimilation",
    "h_b_3"      => 0,
    "k_d_3"      => 0,
    "k_k_3"      => 0,
    "m_n_3"      => 0,
    "s_s_3"      => 0,
    "s_d_3"      => 0,
    "A_eq_3"     => 0,
    "B_eq_3"     => 0,
    "C_eq_3"     => 0
)
    
    
DL1 = Dict(
    "species" => "DL_AMP",
    "cv"      => 0.1,
    
    # AMP parameters
    "E_Hb"       => 0.01746, 
    "E_Hp"       => 1.401,
    "p_M"        => 5686.83,
    "E_G"        => 4400.59,
    "k_J_rate"   => 0.002,
    "v_rate"     => 0.020455,
    "kap"        => 0.331941,
    "kap_R"      => 0.95,
    "shape_factor"       => 0.18946,
    "mortality_constant" => 0.93567,
    "sG"  => 0.001,
    "h_a" => 8.603e-5,
    "L0"  => 0.0001,  

    # food-dependent parameters
    "p_Am" => 651.572, # 2020-03-24 J/d*cm^2
    "y_EX" => 71312.7, # yield of reserve (J) on food (# cells/cm^3)
    "F_m"  => 606.58, # filtration rate (cm^3/day*cm^2)
    
    # embryonal parameters
    "egg_weight" => 3.771781921386719E-4,
    "L_embryo"   => 0.023340451224020977,
    "U_E_embryo" => 6.039673323323045E-5,
    
    # other
    "time_between_molts" => 2.324,
    "crit_mass"          => 0.442002, 
    "starvation_hazard_rate" => 0.784624,
    
    # Toxicity parameters
    # compound 1
    "compound_1" => "Copper", 
    "pmoa_1"     => "Growth", 
    "h_b_1"      => 0/24, # background hazard rate
    "k_d_1"      => 0, # dominant rate constant
    "k_k_1"      => 0, # ?
    "m_n_1"      => 0, # ?
    "s_s_1"      => 0,
    "s_d_1"      => 0, # ?
    "A_eq_1" => 0,
    "B_eq_1" => 0,
    "C_eq_1" => 0,
    
    
    # compound 2
    "compound_2" => "Nickel",
    "pmoa_2"     => "Maintenance",
    "h_b_2"      => 0,
    "k_d_2"      => 0,
    "k_k_2"      => 0,
    "m_n_2"      => 0,
    "s_s_2"      => 0,
    "s_d_2"      => 0,
    "A_eq_2"     => 0,
    "B_eq_2"     => 0,
    "C_eq_2"     => 0,
    
    # compound 3
    "compound_3" => "Zinc",
    "pmoa_3"     => "Assimilation",
    "h_b_3"      => 0,
    "k_d_3"      => 0,
    "k_k_3"      => 0,
    "m_n_3"      => 0,
    "s_s_3"      => 0,
    "s_d_3"      => 0,
    "A_eq_3"     => 0,
    "B_eq_3"     => 0,
    "C_eq_3"     => 0
) 

DL1_ECOREL = Dict(
    "species" => "DL_AMP",
    "cv"      => 0.1,
    
    # AMP parameters
    "E_Hb"       => 0.01746, 
    "E_Hp"       => 1.401,
    "p_M"        => 5686.83,
    "E_G"        => 4400.59,
    "k_J_rate"   => 0.002,
    "v_rate"     => 0.020455,
    "kap"        => 0.331941,
    "kap_R"      => 0.95,
    "shape_factor"       => 0.18946,
    "sG"  => 0.001,
    "h_a" => 8.603e-5,
    "L0"  => 0.0001,  

    # food-dependent parameters
    "p_Am" => 651.572, # 2020-03-24 J/d*cm^2
    "y_EX" => 71312.7, # yield of reserve (J) on food (# cells/cm^3)
    "F_m"  => 606.58, # filtration rate (cm^3/day*cm^2)
    
    # embryonal parameters
    "egg_weight" => 3.771781921386719E-4,
    "L_embryo"   => 0.023340451224020977,
    "U_E_embryo" => 6.039673323323045E-5,
    
    # reproductive timing
    "time_between_molts" => 2.324,
    
    # starvation
    "mortality_constant" => 0.93567,
    "starvation_hazard_rate" => 0.784624,
    "crit_mass" => 0.442002, 
    
    # Toxicity parameters
    # compound 1
    "compound_1" => "Copper", 
    "pmoa_1"     => "Growth", 
    "h_b_1"      => 0/24, # background hazard rate
    "k_d_1"      => 0, # dominant rate constant
    "k_k_1"      => 0, # ?
    "m_n_1"      => 0, # ?
    "s_s_1"      => 0,
    "s_d_1"      => 0, # ?
    "A_eq_1" => 0,
    "B_eq_1" => 0,
    "C_eq_1" => 0,
    
    
    # compound 2
    "compound_2" => "Nickel",
    "pmoa_2"     => "Maintenance",
    "h_b_2"      => 0,
    "k_d_2"      => 0,
    "k_k_2"      => 0,
    "m_n_2"      => 0,
    "s_s_2"      => 0,
    "s_d_2"      => 0,
    "A_eq_2"     => 0,
    "B_eq_2"     => 0,
    "C_eq_2"     => 0,
    
    # compound 3
    # Zinc toxicity parameters based on mean µg dissolved Zn/L
    "compound_3" => "Zinc",
    "pmoa_3"     => "Maintenance",
    "k_d_3"      => 0.727032,
    "m_n_3"      => 344.759,
    "k_k_3"      => 2.67,
    "s_s_3"      => 0.01,
    "A_eq_3"     => 298.85,
    "B_eq_3"     => 2.11,
    "C_eq_3"     => 9.40,
    
    # TODO: remove these parameters entirely from all code
    "h_b_3" => 0,
    "s_d_3" => 0
)