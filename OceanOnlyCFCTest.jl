# Import required packages
using DataFrames, CSV
using DifferentialEquations 
using PythonPlot, Interpolations, NCDatasets
include("source.jl")

datadir(x) = "./" * x
ds = NCDataset(datadir("Tracer_atmospheric_histories_revised_2023.nc"),"r")
ds_CM4X = NCDataset(datadir("Southern_Ocean_Surface_Properties_TimeSeries.nc"),"r")

struct CFCParams
   T
   f
   S::Float64
   P_sfc::Float64
   u::Float64
   compound::String
end

atmospheric_cfc_interp = linear_interpolation(ds["time"][:], ds["CFC11_SH"][:])
 
function atmospheric_cfc(t)
   if t < 2010
      return atmospheric_cfc_interp(t)
   else 
      return atmospheric_cfc_interp(2010)

   end
end


sfc_temperature_forced = linear_interpolation(ds_CM4X["year"][:], ds_CM4X["tos"][1, :])
sfc_temperature_ctrl = linear_interpolation(ds_CM4X["year"][:], ds_CM4X["tos"][2, :])

sfc_si_forced = linear_interpolation(ds_CM4X["year"][:], ds_CM4X["siconc"][1, :])
sfc_si_ctrl = linear_interpolation(ds_CM4X["year"][:], ds_CM4X["siconc"][2, :])

function southern_ocean_temperature(t, scenario)
   if scenario in [:warming, :warming_melt]
      return sfc_temperature_forced(t)
  else
      return sfc_temperature_ctrl(t)
  end
end
function sea_ice_fraction(t, scenario)
   if scenario in [:melt, :warming_melt]
      return sfc_si_forced(t)
   else
      return sfc_si_ctrl(t)
   end
end

# Prevent negative concentrations
function affect!(integrator)
   integrator.u[1] = max(0.0, integrator.u[1])
end

# Check for negative concentrations
function condition(u, t, integrator)
   return u[1] < 0.0 
end

cb = DiscreteCallback(condition, affect!)

# Simplified ODE for ocean only
function cfc_ocean!(du, u, p::CFCParams, t)
   t_real = t + 1950
   CFCocn = u[1] 
   CFCatm = atmospheric_cfc(t_real) * 1e-12
   
   flux = airseagasflux(p.T(t_real), p.S, p.P_sfc, p.u, p.f(t_real), CFCatm, CFCocn) #mol / m² / s
   flux = 3.154e+7  * flux / 5. #mol / m³ / year
   du[1] = flux #mol / m³ / year
end

function run_scenario(scenario, tspan)
   u₀ = [0.0]  # Only ocean initial condition
   
   S = 34.0
   u = 10.0
   
   T = t -> southern_ocean_temperature(t, scenario)
   f = t -> sea_ice_fraction(t, scenario)
   
   p = CFCParams(T, f, S, 1.0, u, "CFC-11")
   prob = ODEProblem(cfc_ocean!, u₀, tspan, p)
   solve(prob, Tsit5(), callback=cb,reltol = 1e-16, abstol = 1e-16)
end

tspan = (0.0, 70.0)
t_plot = range(tspan[1], tspan[2], length= 100)
years = t_plot .+ 1950
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

# Create multi-panel plot
fig, ax = subplots(4, 1, figsize=(10, 10), sharex = true)

# Plot temperature
for scenario in [:warming, :control]
   temp_data =  [southern_ocean_temperature(t, scenario) for t in years]
   ax[0].plot(years, temp_data,
            color=pcolors[scenario], label=plabels[scenario])
   ax[0].annotate(plabels[scenario], 
            xy=(years[end], temp_data[end]),
            xytext=(5, 0), textcoords="offset points",
            va="center", color=pcolors[scenario])
end
ax[0].set_ylabel("Temperature [°C]")
schmidt_number(20)
# Plot sea ice fraction
for scenario in [:melt, :control]
   sea_ice_data = [sea_ice_fraction(t, scenario) for t in years]
   ax[1].plot(years, sea_ice_data,
   color=pcolors[scenario], label=plabels[scenario])
   ax[1].annotate(plabels[scenario], 
                  xy=(years[end], sea_ice_data[end]),
                  xytext=(5, 0), textcoords="offset points",
                  va="center", color=pcolors[scenario])
end
ax[1].set_ylabel("Sea Ice Fraction")
ax[1].set_ylim(0.0, .3)

# Plot atmospheric CFC
ax[2].plot(years, [atmospheric_cfc(t) for t in years],
   color="k", label="Data")
ax[2].set_ylabel("Atmospheric CFC [ppt]")

# Plot ocean CFC with annotations
for scenario in scenarios
   data = ocean_data[scenario] * (1e12 / 1035)
   line = ax[3].plot(years, data, color=pcolors[scenario])
   ax[3].annotate(plabels[scenario], 
                  xy=(years[end], data[end]),
                  xytext=(5, 0), textcoords="offset points",
                  va="center", color=pcolors[scenario])
end

ax[3].set_ylabel("Ocean CFC [pmol/m³]")

# Add grids and legends to all panels
[a.grid(alpha = 0.4) for a in ax]
[a.legend(ncols = 2) for a in ax]

# Set x-axis properties
ax[3].set_xlabel("Year")
ax[3].set_xlim(1950, 2100)
ax[3].set_xticks(iyears[1:20:end])
fig.tight_layout()
fig