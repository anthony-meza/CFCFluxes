using DataFrames, CSV, 
PythonPlot, Unitful
using UnitfulLinearAlgebra
ENV["UNITFUL_FANCY_EXPONENTS"] = true

include("../src/src.jl")

function meshgrid(x, y)
    return (x' .* ones(length(y)), ones(length(x))' .* y)
end

t = -1:0.1:40 
s = 0:0.1:40
T, S = meshgrid(t, s) .* (u"°C", u"g/kg")

F = calculate_solubility.(T, S; compound = "CFC-11")

fig, ax = subplots()
cb = ax.contourf(ustrip(T), ustrip(S), ustrip(F), cmap = "Spectral_r", levels = 20)
fig.colorbar(cb, label = "mol atm⁻¹ m⁻³")
ax.set_title("Solubility of CFC-11 according to\n Orr et. al. (2017)")
ax.set_xlabel("Temperature [deg C]")
ax.set_ylabel("Salinity [g/kg]")
fig