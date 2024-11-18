#= 
Core CFC functions implementing gas exchange calculations following:
Orr, J. et al. (2017). 
Biogeochemical protocols and diagnostics for the CMIP6 Ocean Model 
Intercomparison Project (OMIP). Geoscientific Model Development, 10(6), 2169–2199.
https://doi.org/10.5194/gmd-10-2169-2017
=# 

using DataFrames, CSV, Unitful, UnitfulMoles

#Defining Units#
ppm = u"µmol/mol" #number molucules parts per million number molecules of dry air
ppb = u"nmol/mol" #number molucules parts per billion number molecules of dry air
ppt = u"pmol/mol" #number molucules parts per trillion number molecules of dry air
degC = u"°C"
molm3 = u"mol/m^3"
rho = 1024u"kg/m^3"

#Defining Approximation Coefficients#
solubility_types = Dict(:Coefficients => String, :CFC11 => Float32, :CFC12 => Float32)
solubility_coefficients = DataFrame(CSV.File(datadir("SolubilityCoefficients.csv"); header=2, types = solubility_types))

schmidt_number_types = Dict(:Coefficients => String, :CFC11 => Float32, :CFC12 => Float32)
schmidt_number_coefs = DataFrame(CSV.File(datadir("SchmidtNumber.csv"); header=2, types = schmidt_number_types))

#Defining Unit Checks#
check_TdegC_units(T) = (unit(T) != degC) && error("T is not in units Celsius")
check_S_units(S) = (unit(S) != u"g/kg") && error("S is not in units grams per kilogram of seawater")
check_u_units(u) = (unit(u) != u"m/s") && error("u is not in units meters per second")
check_P_units(P) = (unit(P) != u"atm") && error("P is not in units atmoshpheres")


"""
   calculate_solubility(T::Unitful.Temperature, S::Any; compound::String="CFC-11")

Calculate the solubility coefficient (F) of a CFC compound in seawater.

# Arguments
- `T`: Temperature (must be in degrees Celsius)
- `S`: Salinity (g/kg)
- `compound`: ("CFC-11" or "CFC-12")

# Returns
- Solubility coefficient in mol/(m³⋅atm)
"""
function calculate_solubility(T::Unitful.Temperature, 
                             S::Any; 
                             compound::String="CFC-11")

    check_TdegC_units(T); check_S_units(S)

    T_kelvin = uconvert(u"K", T)
    coefficient_column = compound == "CFC-11" ? 2 : 
                        compound == "CFC-12" ? 3 : 
                        error("Unsupported compound: $compound")
    
    coefficients = solubility_coefficients[:, coefficient_column]
    
    α₁, α₂, α₃, α₄ = coefficients[1:4]
    β₁, β₂, β₃ = (coefficients[5:7] .* u"kg/g")

    T_ratio = T_kelvin / 100u"K"
    α_terms = α₁ + (α₂ * inv(T_ratio)) + (α₃ * log(T_ratio)) + (α₄ * (T_ratio^2))
    
    β_terms = S * (β₁ + (β₂ * T_ratio) + (β₃ * (T_ratio^2)))
    lnF = α_terms + β_terms
    expF = exp(lnF) * 1u"mol/L/atm"

    F = uconvert(u"mol/((m^3) * atm)", expF)
    return F
end

"""
   schmidt_number(T::Unitful.Temperature; compound::String="CFC-11")

Calculate the Schmidt number for a CFC compound.

# Arguments
- `T`: Temperature (must be in °C)
- `compound`: ("CFC-11" or "CFC-12")

# Returns
- Dimensionless Schmidt number
"""
function schmidt_number(T::Unitful.Temperature; compound::String="CFC-11")

    check_TdegC_units(T)

    coefficient_column = compound == "CFC-11" ? 2 : 
                        compound == "CFC-12" ? 3 : 
                        error("Unsupported compound: $compound")
    
    coefficients = schmidt_number_coefs[:, coefficient_column]
    A, B, C, D, E = coefficients[1:5] .* [u"K"^(-j) for j in 0:4]

    T_a = T - 0degC #convert to temperature difference to get this to work with Unitful

    Sc = A + (B * T_a) + (C * T_a^2) + (D * T_a^3) + (E * T_a^4)
    return Sc
end

"""
   piston_velocity(T::Unitful.Temperature, u::Unitful.Velocity, f; compound::String="CFC-11")

Calculate the air-sea gas exchange piston velocity.

# Arguments
- `T`: Temperature (must be in °C)
- `u`: Wind speed (m/s)
- `f`: Ice fraction (0-1)
- `compound`: ("CFC-11" or "CFC-12")

# Returns
- Piston velocity
"""
function piston_velocity(T::Unitful.Temperature, 
                        u::Unitful.Velocity, f; compound::String="CFC-11")
    check_TdegC_units(T); check_u_units(u)

    a = 6.97e-7u"s/m" 
    Sc = schmidt_number(T; compound) #no units
    ScC02 = 660 #unitless
    kw = a * inv(sqrt(Sc / ScC02)) * (u^2) * (1 - f)
    return kw
end


"""
   gassat(T, S, P_sfc, x_compound; compound::String="CFC-11")

Calculate the saturation concentration of a CFC in seawater.

# Arguments
- `T`: Temperature (must be in °C)
- `S`: Salinity (g/kg)
- `P_sfc`: Surface pressure (atm)
- `x_compound`: CFC mixing ratio in ppt
- `compound`: CFC compound type ("CFC-11" or "CFC-12")

# Returns
- Saturation concentration in mol/m^3
"""
function gassat(T, S, P_sfc, x_compound; compound::String="CFC-11")
    check_TdegC_units(T); check_S_units(S); check_P_units(P_sfc); 
    (unit(x_compound) != ppt) && error("P is not in units parts per trillion")

    P₀ = 1u"atm" #reference atmosphere 
    F = calculate_solubility(T, S; compound)

    #convert from a mixing ratio (i.e. ppt, mass fraction) to an atmospheric partial pressure
    p_compound = uconvert(u"mol/mol", x_compound)u"atm"  #percentage of atmospheric pressure from compound
    # p_compound = uconvert(u"patm", p_compound)
    Xsat = P_sfc * F * (p_compound) / P₀
    return Xsat
end

"""
   gasatmconc(T, S, P_sfc, Xsat; compound::String="CFC-11")

Convert dissolved CFC concentration to atmospheric mixing ratio.

# Arguments
- `T`: Temperature (must be in °C)
- `S`: Salinity (g/kg)
- `P_sfc`: Surface pressure (atm)
- `Xsat`: Dissolved CFC concentration in mol/m³
- `compound`: ("CFC-11" or "CFC-12")

# Returns
- Atmospheric mixing ratio in ppt
"""
function gasatmconc(T, S, P_sfc, Xsat; compound::String="CFC-11")
    check_TdegC_units(T); check_S_units(S); check_P_units(P_sfc); 
    (unit(Xsat) != molm3) && error("P is not in units parts per trillion")

    P₀ = 1u"atm" #reference atmosphere 
    F = calculate_solubility(T, S; compound)
    p_compound = (Xsat * P₀) / (P_sfc * F)

    #convert from a mixing ratio (i.e. ppt, mass fraction) to an atmospheric partial pressure
    x_compound = uconvert(u"pmol/mol", ustrip(p_compound)) #percentage of atmospheric pressure from compound
    return x_compound
end

"""
   airseagasflux(T, S, P_sfc, u, f, xcfc, CFCocn; compound::String="CFC-11")

Calculate the air-sea gas flux of CFCs.

# Arguments
- `T`: Temperature (must be in °C)
- `S`: Salinity  
- `P_sfc`: Surface pressure
- `u`: Wind speed
- `f`: Ice fraction (0-1)
- `xcfc`: Atmospheric CFC concentration in ppt
- `CFCocn`: Oceanic CFC concentration in mol/m^3
- `compound`: CFC compound type ("CFC-11" or "CFC-12")

# Returns
- Air-sea gas flux (mol/m^3/s)
"""
function airseagasflux(T, S, P_sfc, u, f, xcfc, CFCocn; compound::String="CFC-11")
    check_TdegC_units(T); check_S_units(S); check_P_units(P_sfc); 
    (unit(xcfc) != ppt) && error("xcfc is not in units parts per trillion, x is:" * string(xcfc))

    k = piston_velocity(T, u, f; compound)


    Asat = gassat(T, S, P_sfc, xcfc; compound)
    return k * (Asat - CFCocn)
end