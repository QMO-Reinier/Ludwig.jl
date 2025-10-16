using PlotlyJS
using LaTeXStrings
using LinearAlgebra
include("Sr2RuO4.jl")

E_F = 0

# We make arrays of k values, from the center of the BZ to the edge
k_array_x = range(-0.5,0.5,201)
k_array_y = range(-0.5,0.5,201)

# Create empty arrays of correct size for all bands
band_alpha = zeros(Float64,201,201)
band_beta = zeros(Float64,201,201)
band_gamma = zeros(Float64,201,201)

# This loop uses functions from Sr2RuO4.jl to calculate dispersion relations
for (i, k_x) in enumerate(k_array_x)
    for (j, k_y) in enumerate(k_array_y)
        band_alpha[i,j] = ham_α([k_x,k_y])
        band_beta[i,j] = ham_β([k_x,k_y])
        band_gamma[i,j] = ham_γ([k_x,k_y])
    end
end

# Plot contours
ticks_x      = [-0.5,0,0.5]
ticklabels_x = ["-π/a",0,"π/a"]
ticks_y      = [0,0.5]
ticklabels_y = [0,"π/a"]

traces = [
    contour(
        x = k_array_x, y = k_array_y, z = band_alpha,
        autocontour=false,
        contours_start=E_F, contours_end=E_F, contours_size=1e-5,
        contours_coloring="lines",
        name="α", line=attr(width=3),
        colorscale=[[0,"blue"], [1,"blue"]],
        showlegend=true, showscale=false
    ),
    contour(
        x = k_array_x, y = k_array_y, z = band_beta,
        autocontour=false,
        contours_start=E_F, contours_end=E_F, contours_size=1e-5,
        contours_coloring="lines",
        name="β", line=attr(width=3),
        colorscale=[[0,"green"], [1,"green"]],
        showlegend=true ,showscale=false
    ),
    contour(
        x = k_array_x, y = k_array_y, z = band_gamma,
        autocontour=false,
        contours_start=E_F, contours_end=E_F, contours_size=1e-5,
        contours_coloring="lines",
        name="γ", line=attr(width=3),
        colorscale=[[0,"red"], [1,"red"]],
        showlegend=true, showscale=false
    )
]

p = Plot(traces, Layout(
    xaxis=attr(title="k_x", tickvals=ticks_x, ticktext=ticklabels_x, scaleanchor="y", constrain="domain"),
    yaxis=attr(title="k_y", tickvals=ticks_y, ticktext=ticklabels_y),
    width=600, height=600,
    font=attr(size=20),
))

p