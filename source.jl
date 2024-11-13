using DataFrames, CSV

"""
    calculate_solubility(T, S; compound = "CFC-11")

Calculate the solubility coefficient (F) for specified CFC compounds in seawater.

# Arguments
- `T::Number`: Temperature in degrees Celsius
- `S::Number`: Salinity in practical salinity units (PSU)
- `compound::String`: Chemical compound, either "CFC-11" or "CFC-12"

# Returns
- `Float64`: Solubility coefficient F in mol/kg/atm

# Details
Calculates the solubility coefficient F using the equation:
ln(F) = a₁ + (a₂ * 100/T) + (a₃ * ln(T/100)) + (a₄ * (T/100)²) + 
        S * (a₅ + (a₆ * T/100) + (±a₇ * (T/100)²))

where T is in Kelvin and a₁-a₇ are empirically determined coefficients.
Note: The a₇ term is negative for CFC-12 and positive for CFC-11.

# References
Warner, M. J., & Weiss, R. F. (1985). Solubilities of chlorofluorocarbons 11 and 12 in
water and seawater. Deep Sea Research Part A, 32(12), 1485-1497.

# Examples
```julia
F = calculate_solubility(25.0, 35.0, compound="CFC-11")
F = calculate_solubility(20.0, 34.5, compound="CFC-12")
```
"""
function calculate_solubility(T, S; compound::String = "CFC-11")
    # Convert temperature from Celsius to Kelvin
    T_kelvin = T + 273.15
    
    # Read coefficients from file
    coefficient_column = compound == "CFC-11" ? 2 : 
                        compound == "CFC-12" ? 3 : 
                        error("Unsupported compound: $compound")
    
    coefficients = DataFrame(CSV.File("SolubilityCoefficients.csv"; header=2))[:, coefficient_column]
    
    # Calculate intermediate terms for readability
    T_ratio = T_kelvin / 100
    temp_terms = coefficients[1] + 
                (coefficients[2] * 100 / T_kelvin) + 
                (coefficients[3] * log(T_ratio)) + 
                (coefficients[4] * T_ratio^2)
    
    # Note: CFC-12 uses negative a₇, CFC-11 uses positive a₇
    salinity_terms = S * (coefficients[5] + 
                         (coefficients[6] * T_ratio) + 
                         (coefficients[7] * T_ratio^2))
    
    # Calculate ln(F) and return F
    lnF = temp_terms + salinity_terms
    return exp(lnF)
end