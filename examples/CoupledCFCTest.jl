# Import required packages
using DataFrames, CSV, DifferentialEquations 
using PythonPlot, Interpolations, NCDatasets
include("../src/src.jl")

# Load data
ds = NCDataset(datadir("Tracer_atmospheric_histories_revised_2023.nc"), "r")
ds_CM4X = NCDataset(datadir("Southern_Ocean_Surface_Properties_TimeSeries.nc"), "r")
ds_CM4X
# Create interpolations
atmospheric_cfc_interp = linear_interpolation(ds["time"][:], ds["CFC11_SH"][:])
sfc_temperature_forced = linear_interpolation(ds_CM4X["year"][:], ds_CM4X["tos"][1, :])
sfc_temperature_ctrl = linear_interpolation(ds_CM4X["year"][:], ds_CM4X["tos"][2, :])
sfc_si_forced = linear_interpolation(ds_CM4X["year"][:], ds_CM4X["siconc"][1, :])
sfc_si_ctrl = linear_interpolation(ds_CM4X["year"][:], ds_CM4X["siconc"][2, :])

function southern_ocean_temperature(t, scenario)
   scenario in [:warming, :warming_melt] ? sfc_temperature_forced(t) : sfc_temperature_ctrl(t)
end

function sea_ice_fraction(t, scenario)
  scenario in [:melt, :warming_melt] ? sfc_si_forced(t) : sfc_si_ctrl(t)
end

# Define parameters struct for CFC calculations
struct CFCParams
   T
   f
   S
   P_sfc
   u
   compound::String
end

# System of ODEs for coupled ocean-atmosphere CFC evolution

function cfc_coupled_system!(dx, x, p::CFCParams, t)
   t_real = t + 1950

   CFCocn, CFCatmsat = x .* molm3

   CFCatm = gasatmconc(p.T(t_real), p.S, p.P_sfc,CFCatmsat)
   flux = airseagasflux(p.T(t_real), p.S, p.P_sfc, p.u, p.f(t_real), CFCatm, CFCocn) #mol / m² / s

   dz = 5u"m" #depth of surface layer
   uflux =  uconvert(u"mol/m^3/yr", flux / dz) #mol / m³ / year
   dx[1] = ustrip(uflux)
   anth_flux = gassat(p.T(t_real), p.S, p.P_sfc, anthropogenic_source(t)ppt)
   dx[2] = -ustrip(uflux) + ustrip(anth_flux)

end


# Run simulation for a given scenario
function run_scenario(scenario, tspan)
   S = 34.0u"g/kg"
   uvel = 15.0u"m/s"
   Psfc = 1.0u"atm"
   
   p = CFCParams(
       t -> southern_ocean_temperature(t, scenario)degC,
       t -> sea_ice_fraction(t, scenario),
       S, Psfc, uvel, "CFC-11"
   )

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


   u₀ = [0.0, 0.0]    # Initial conditions
   prob = ODEProblem(cfc_coupled_system!, u₀, tspan, p)
   solve(prob, Tsit5(), callback=cb,reltol = 1e-16, abstol = 1e-16)
end

# Run simulations
tspan = (0.0, 140.0)
t_plot = range(tspan[1], tspan[2], length= 100)
years = range(tspan[1], tspan[2], length=100) .+ 1950
iyears = years[1]:1:years[end]
scenarios = [:control, :warming, :melt, :warming_melt]
solutions = Dict(scenario => run_scenario(scenario, tspan) for scenario in scenarios)

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
for scenario in [:warming, :control]
   ax[0].plot(years, [southern_ocean_temperature(t, scenario) for t in years],
            color=pcolors[scenario], label=plabels[scenario])
end
ax[0].set_ylabel("Temperature [°C]")

# Plot sea ice fraction
for scenario in [:melt, :control]
   ax[1].plot(years, [sea_ice_fraction(t, scenario) for t in years],
   color=pcolors[scenario], label=plabels[scenario])
end
ax[1].set_ylabel("Sea Ice Fraction")

# Plot atmospheric CFC
for scenario in scenarios
   ax[2].plot(years, atm_data[scenario],
   color=pcolors[scenario], label=plabels[scenario])
end
ax[2].set_ylabel("Atmospheric CFC Equillibrium concentration [mol/m^3]")
fig
# Plot ocean CFC with annotations
for scenario in scenarios
   line = ax[3].plot(years, ocean_data[scenario], color=pcolors[scenario])
   ax[3].annotate(plabels[scenario], 
                  xy=(years[end], ocean_data[scenario][end]),
                  xytext=(5, 0), textcoords="offset points",
                  va="center", color=pcolors[scenario])
end
ax[3].set_ylabel("Ocean CFC  [mol/m^3]")

# Add grids and legends to all panels
[a.grid(alpha = 0.4) for a in ax]
[a.legend(ncols = 2) for a in ax]

# Set x-axis properties
ax[3].set_xlabel("Year")
ax[3].set_xlim(1950, 2130)
ax[3].set_xticks(iyears[1:20:end])
fig.tight_layout()
fig