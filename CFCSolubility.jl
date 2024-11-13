using DataFrames, CSV, PythonPlot

include("source.jl")

function meshgrid(x, y)
    return (x' .* ones(length(y)), ones(length(x))' .* y)
end

t = -1:0.1:40
s = 0:0.1:40
T, S = meshgrid(t, s)
F = calculate_solubility.(T, S; compound = "CFC-11")

fig, ax = subplots()
cb = ax.contourf(T, S, 100 .* F, cmap = "Spectral_r", levels = 20)
fig.colorbar(cb)
ax.set_title("Solubility of CFC-11 according to\nWarner and Weiss (1985)")
ax.set_xlabel("Temperature [deg C]")
ax.set_ylabel("Salinity [%]")
fig