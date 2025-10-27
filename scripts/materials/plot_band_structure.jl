using PlotlyJS
using LaTeXStrings
using LinearAlgebra
include("Sr2RuO4.jl")

# We make arrays of k values, from the center of the BZ to the edge
k_array_x = range(0,0.5,101)
k_array_y = range(0,0.5,101)

# Create empty arrays of correct size for all bands
band_alpha = zeros(Float64,101,101)
band_beta = zeros(Float64,101,101)
band_gamma = zeros(Float64,101,101)

# This loop uses functions from Sr2RuO4.jl to calculate dispersion relations
for (i, k_x) in enumerate(k_array_x)
    for (j, k_y) in enumerate(k_array_y)
        band_alpha[i,j] = ham_α([k_x,k_y])
        band_beta[i,j] = ham_β([k_x,k_y])
        band_gamma[i,j] = ham_γ([k_x,k_y])
    end
end

# Take paths across BZ to vizualize bandstructure
band_alpha_gamma_to_X = band_alpha[:,1]
band_alpha_X_to_M     = band_alpha[end,:]
band_alpha_M_to_gamma = reverse(diag(band_alpha))
band_beta_gamma_to_X  = band_beta[:,1]
band_beta_X_to_M      = band_beta[end,:]
band_beta_M_to_gamma  = reverse(diag(band_beta))
band_gamma_gamma_to_X = band_gamma[:,1]
band_gamma_X_to_M     = band_gamma[end,:]
band_gamma_M_to_gamma = reverse(diag(band_gamma))

# Concatenate the three paths
total_band_alpha = vcat(band_alpha_gamma_to_X,band_alpha_X_to_M,band_alpha_M_to_gamma)
total_band_beta  = vcat(band_beta_gamma_to_X,band_beta_X_to_M,band_beta_M_to_gamma)
total_band_gamma = vcat(band_gamma_gamma_to_X,band_gamma_X_to_M,band_gamma_M_to_gamma)

# Plotting bandstructure
ticks = [0, length(band_alpha_gamma_to_X), length(band_alpha_gamma_to_X)+length(band_alpha_X_to_M),
        length(band_alpha_gamma_to_X)+length(band_alpha_X_to_M)+length(band_alpha_M_to_gamma)]
ticklabels = ["Γ", "X", "M", "Γ"]

plt =    plot(
            [
                scatter(y=total_band_alpha, mode="lines", name="α band", line=attr(color="blue")),
                scatter(y=total_band_beta,  mode="lines", name="β band", line=attr(color="green")),
                scatter(y=total_band_gamma, mode="lines", name="γ band", line=attr(color="red"))
            ],
            Layout(
                xaxis=attr(
                    tickvals=ticks,
                    ticktext=ticklabels
                ),
                yaxis_title="E(k)",
                font=attr(size=20),
                legend=attr(title="Bands"),
            )
        )

display(plt)

bands = [band_alpha,band_beta,band_gamma]