using DataFrames, CSV

# ===== Core CFC functions ===== #
function calculate_solubility(T, S; compound::String="CFC-11")
    types = Dict(
    :Coefficients => String,
    :CFC11 => Float32, 
    :CFC12 => Float32
    )

    T_kelvin = T + 273.15
    coefficient_column = compound == "CFC-11" ? 2 : 
                        compound == "CFC-12" ? 3 : 
                        error("Unsupported compound: $compound")
    
    coefficients = DataFrame(CSV.File("SolubilityCoefficients.csv"; header=2, types = types))
    coefficients = coefficients[:, coefficient_column]
    
    T_ratio = T_kelvin / 100
    temp_terms = coefficients[1] + 
                (coefficients[2] * 100 / T_kelvin) + 
                (coefficients[3] * log(T_ratio)) + 
                (coefficients[4] * T_ratio^2)
    
    salinity_terms = S * (coefficients[5] + 
                         (coefficients[6] * T_ratio) + 
                         (coefficients[7] * T_ratio^2))
    
    return exp(temp_terms + salinity_terms)
end

function schmidt_number(T; compound::String="CFC-11")
    types = Dict(
    :Coefficients => String,
    :CFC11 => Float32, 
    :CFC12 => Float32
    )


    coefficient_column = compound == "CFC-11" ? 2 : 
                        compound == "CFC-12" ? 3 : 
                        error("Unsupported compound: $compound")
    
    coefficients = DataFrame(CSV.File("SchmidtNumber.csv"; header=2, types = types))
    coefficients = coefficients[:, coefficient_column]

    return coefficients[1] + 
           (coefficients[2] * T) + 
           (coefficients[3] * T^2) + 
           (coefficients[4] * T^3) + 
           (coefficients[5] * T^4)
end

function piston_velocity(T, u, f; compound::String="CFC-11")
    a = 6.97e-7
    Sc = schmidt_number(T; compound)
    return a * inv(sqrt(Sc / 660)) * u^2 * (1 - f)
end

function gassat(T, S, P_sfc, xcfc; compound::String="CFC-11")
    P₀ = 1    
    F = calculate_solubility(T, S; compound)
    return P_sfc * F * (xcfc * 1e-12) / P₀
end

function airseagasflux(T, S, P_sfc, u, f, xcfc, CFCocn; compound::String="CFC-11")
    k = piston_velocity(T, u, f; compound)
    Asat = gassat(T, S, P_sfc, xcfc; compound)
    return k * (Asat - CFCocn)
end

# Scenario functions
function southern_ocean_temperature(t, scenario)
    year = t + 1950
    base_temp = 0.0
    seasonal_amp = 0.0
    seasonal = seasonal_amp * sin(2π * (t % 1))
    
    if scenario in [:warming, :warming_melt]
        warming_rate = 0.03  # °C/year
        warming = warming_rate * (year - 1950)
        return base_temp + warming + seasonal
    else
        return base_temp + seasonal
    end
end

function sea_ice_fraction(t, scenario)
    year = t + 1950
    base_ice = 0.5
    seasonal_amp = 0.0
    seasonal = seasonal_amp * sin(2π * (t % 1))
    
    if scenario in [:melt, :warming_melt]
        melt_rate = 0.004  # per year
        ice = base_ice - melt_rate * (year - 1950) + seasonal
        return clamp(ice, 0.0, 1.0)
    else
        ice = base_ice + seasonal
        return clamp(ice, 0.0, 1.0)
    end
end

function anthropogenic_source(t)
    year = t + 1950
    
    if year < 1995
        return 10.0 * (1 - exp(-0.1 * (year - 1950)))
    elseif year < 2010
        decay_rate = 0.1
        return 15.0 * exp(-decay_rate * (year - 1995))
    else
        return 0.0
    end
end