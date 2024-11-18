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

# Time-dependent functions
function atmospheric_cfc(t)
    t < 2010 ? atmospheric_cfc_interp(t) : atmospheric_cfc_interp(2010)
end

function southern_ocean_temperature(t, scenario)
    scenario in [:warming, :warming_melt] ? sfc_temperature_forced(t) : sfc_temperature_ctrl(t)
end

function sea_ice_fraction(t, scenario)
   scenario in [:melt, :warming_melt] ? sfc_si_forced(t) : sfc_si_ctrl(t)
end

# Parameter struct
struct CFCParams
   T
   f
   S
   P_sfc
   u
   compound::String
end

# Simplified ODE for ocean only
function cfc_ocean!(dx, x, p::CFCParams, t)
   t_real = t + 1950

   CFCocn = x[1]molm3
   CFCatm = atmospheric_cfc(t_real)ppt
   
   flux = airseagasflux(p.T(t_real), p.S, 
                        p.P_sfc, p.u, 
                        p.f(t_real), 
                        CFCatm, CFCocn) #mol / m² / s
   dz = 5u"m" #depth of surface layer
   uflux =  uconvert(u"mol/m^3/yr", flux / dz) #mol / m³ / year
   dx[1] = ustrip(uflux)
end

function run_scenario(scenario, tspan)

   S = 34.0u"g/kg"
   uvel = 10.0u"m/s"
   Psfc = 1.0u"atm"
   
   p = CFCParams(
       t -> southern_ocean_temperature(t, scenario)degC,
       t -> sea_ice_fraction(t, scenario),
       S, Psfc, uvel, "CFC-11"
   )

   # Prevent negative concentrations
   condition(u, t, integrator) = u[1] < 0.0 
   affect!(integrator) = (integrator.u[1] = max(0.0, integrator.u[1]))
   cb = DiscreteCallback(condition, affect!)
   u₀ = [0.0]  # Only ocean initial condition
   prob = ODEProblem(cfc_ocean!, u₀, tspan, p)
   solve(prob, Tsit5(), callback=cb,reltol = 1e-16, abstol = 1e-16)
end


# Run simulations
tspan = (0.0, 70.0)
t_plot = range(tspan[1], tspan[2], length= 100)
years = range(tspan[1], tspan[2], length=100) .+ 1950
iyears = years[1]:1:years[end]
scenarios = [:control, :warming, :melt, :warming_melt]
solutions = Dict(scenario => run_scenario(scenario, tspan) for scenario in scenarios)

# Plotting setup

colors = Dict(:control => "k", :warming => "r", :melt => "b", :warming_melt => "purple")
labels = Dict(:control => "Control", :warming => "Warming", :melt => "Ice Melt", :warming_melt => "Warming + Melt")

# Process data for plotting
ocean_data = Dict(
    scenario => uconvert.(
        u"pmol/kg", 
        [u[1] for u in sol.(t_plot)]u"mol/m^3" / 1024u"kg/m^3"
    ) for (scenario, sol) in solutions
)

# Create plots
fig, ax = subplots(4, 1, figsize=(10, 10), sharex=true)

# Temperature panel
for scenario in [:warming, :control]
    temp_data = [southern_ocean_temperature(t, scenario) for t in years]
    ax[0].plot(years, temp_data, color=colors[scenario], label=labels[scenario])
    ax[0].annotate(labels[scenario], xy=(years[end], temp_data[end]),
                  xytext=(5, 0), textcoords="offset points",
                  va="center", color=colors[scenario])
end
ax[0].set_ylabel("Temperature [°C]")

# Sea ice panel
for scenario in [:melt, :control]
    ice_data = [sea_ice_fraction(t, scenario) for t in years]
    ax[1].plot(years, ice_data, color=colors[scenario], label=labels[scenario])
    ax[1].annotate(labels[scenario], xy=(years[end], ice_data[end]),
                  xytext=(5, 0), textcoords="offset points",
                  va="center", color=colors[scenario])
end
ax[1].set_ylabel("Sea Ice Fraction")
ax[1].set_ylim(0.0, 0.3)

# Atmospheric CFC panel
ax[2].plot(years, [atmospheric_cfc(t) for t in years], color="k", label="Data")
ax[2].set_ylabel("Atmospheric CFC [ppt]")

# Ocean CFC panel
for scenario in scenarios
   units = unit(first(ocean_data[:warming]))
    data = ustrip(ocean_data[scenario])
    ax[3].plot(years, data, color=colors[scenario])
    ax[3].annotate(labels[scenario], xy=(years[end], data[end]),
                  xytext=(5, 0), textcoords="offset points",
                  va="center", color=colors[scenario])
   ax[3].set_ylabel("Ocean CFC [$units]")

end
# Format plots
[a.grid(alpha=0.4) for a in ax]
[a.legend(ncols=2) for a in ax]
ax[3].set_xlabel("Year")
ax[3].set_xlim(1950, 2100)
ax[3].set_xticks(iyears[1:20:end])
fig.tight_layout()
fig