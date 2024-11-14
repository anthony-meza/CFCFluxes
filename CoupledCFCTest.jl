# Import required packages
using DataFrames, CSV
using DifferentialEquations 
using PythonPlot
include("source.jl")

# Define parameters struct for CFC calculations
struct CFCParams
   T                   # Temperature function
   f                   # Sea ice fraction function
   S::Float64         # Salinity
   P_sfc::Float64     # Surface pressure
   u::Float64         # Wind speed
   compound::String   # CFC compound name
end

# System of ODEs for coupled ocean-atmosphere CFC evolution
function cfc_coupled_system!(du, u, p::CFCParams, t)
   CFCocn, CFCatm = u
   # Calculate air-sea gas exchange
   flux = airseagasflux(p.T(t), p.S, p.P_sfc, p.u, p.f(t), CFCatm, CFCocn)
   # Evolution equations for ocean and atmosphere concentrations
   du[1] = flux / 1024.
   du[2] = (-flux + anthropogenic_source(t)) / 1.
end

# Prevent negative concentrations
function affect!(integrator)
   integrator.u[1] = max(0.0, integrator.u[1])
   integrator.u[2] = max(0.0, integrator.u[2])
end

# Check for negative concentrations
function condition(u, t, integrator)
   return u[1] < 0.0 || u[2] < 0.0
end

cb = DiscreteCallback(condition, affect!)

# Run simulation for a given scenario
function run_scenario(scenario)
   tspan = (0.0, 150.0)
   u₀ = [0.0, 0.0]    # Initial conditions
   
   # Fixed parameters
   S = 34.0           # Salinity
   u = 10.0          # Wind speed
   
   # Time-dependent temperature and ice fraction
   T = t -> southern_ocean_temperature(t, scenario)
   f = t -> sea_ice_fraction(t, scenario)
   
   # Setup and solve ODE problem
   p = CFCParams(T, f, S, 1.0, u, "CFC-11")
   prob = ODEProblem(cfc_coupled_system!, u₀, tspan, p)
   solve(prob, Tsit5(), callback=cb)
end

# Run all scenarios
scenarios = [:control, :warming, :melt, :warming_melt]
solutions = Dict(scenario => run_scenario(scenario) for scenario in scenarios)

# Setup time arrays for plotting
iyears = 1950:2100
years = range(1950, 2100, length=2 * length(iyears))
t_plot = range(0, 150, length= 2 * length(iyears))

# Define plotting colors and labels
pcolors = Dict(:control => "k", :warming => "r", 
            :melt => "b", :warming_melt => "purple")
plabels = Dict(:control => "Control", :warming => "Warming", 
            :melt => "Ice Melt", :warming_melt => "Warming + Melt")

# Extract solution data for plotting
ocean_data = Dict(scenario => [u[1] for u in sol.(t_plot)] for (scenario, sol) in solutions)
atm_data = Dict(scenario => [u[2] for u in sol.(t_plot)] for (scenario, sol) in solutions)

# Create multi-panel plot
fig, ax = subplots(4, 1, figsize=(10, 10), sharex = true)

# Plot temperature
for scenario in scenarios
   ax[0].plot(years, [southern_ocean_temperature(t, scenario) for t in t_plot],
            color=pcolors[scenario], label=plabels[scenario])
end
ax[0].set_ylabel("Temperature [°C]")

# Plot sea ice fraction
for scenario in scenarios
   ax[1].plot(years, [sea_ice_fraction(t, scenario) for t in t_plot],
   color=pcolors[scenario], label=plabels[scenario])
end
ax[1].set_ylabel("Sea Ice Fraction")
ax[1].set_ylim(-0.1, 1.1)

# Plot atmospheric CFC
for scenario in scenarios
   ax[2].plot(years, atm_data[scenario],
   color=pcolors[scenario], label=plabels[scenario])
end
ax[2].set_ylabel("Atmospheric CFC [ppt]")

# Plot ocean CFC with annotations
for scenario in scenarios
   line = ax[3].plot(years, ocean_data[scenario], color=pcolors[scenario])
   ax[3].annotate(plabels[scenario], 
                  xy=(years[end], ocean_data[scenario][end]),
                  xytext=(5, 0), textcoords="offset points",
                  va="center", color=pcolors[scenario])
end
ax[3].set_ylabel("Ocean CFC [mol/kg]")

# Add grids and legends to all panels
[a.grid(alpha = 0.4) for a in ax]
[a.legend(ncols = 2) for a in ax]

# Set x-axis properties
ax[3].set_xlabel("Year")
ax[3].set_xlim(1950, 2130)
ax[3].set_xticks(iyears[1:20:end])
fig.tight_layout()
fig