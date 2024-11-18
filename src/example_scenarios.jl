# Scenario functions

"""
    southern_ocean_temperature(t, scenario)

Calculate the Southern Ocean temperature for a given time point and scenario.

# Arguments
- `t::Float64`: Time in years since 1950
- `scenario::Symbol`: Climate scenario to simulate. Valid options are:
    - `:warming`: Temperature increases linearly
    - `:warming_melt`: Temperature increases with concurrent ice melt
    - other symbols: No temperature change

# Returns
- `Float64`: Temperature in degrees Celsius relative to base temperature

# Details
The temperature evolution consists of:
- Base temperature (0.0°C)
- Seasonal variation (currently disabled, amplitude = 0.0)
- Linear warming of 0.03°C/year since 1950 (only for warming scenarios)
"""
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

"""
    sea_ice_fraction(t, scenario)

Calculate the sea ice fraction in the Southern Ocean for a given time point and scenario.

# Arguments
- `t::Float64`: Time in years since 1950
- `scenario::Symbol`: Climate scenario to simulate. Valid options are:
    - `:melt`: Ice melts linearly
    - `:warming_melt`: Ice melts with concurrent warming
    - other symbols: Constant ice coverage

# Returns
- `Float64`: Sea ice fraction between 0.0 and 1.0

# Details
The ice fraction evolution consists of:
- Base ice fraction (0.5 or 50%)
- Seasonal variation (currently disabled, amplitude = 0.0)
- Linear melting of 0.4%/year since 1950 (only for melting scenarios)
- Output is clamped between 0.0 and 1.0
```
"""
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

"""
    anthropogenic_source(t)

Calculate anthropogenic source contribution for a given time point.
Source is in ppt/year
# Arguments
- `t::Float64`: Time in years since 1950

# Returns
- `Float64`: Source contribution in arbitrary units

# Details
The source evolution has three distinct periods:
1. Growth phase (1950-1994):
   - Exponential approach to 10.0 with rate constant 0.1
2. Decay phase (1995-2009):
   - Exponential decay from 15.0 with rate constant 0.1
3. Zero contribution (2010 onwards):
   - Returns 0.0
```
"""
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