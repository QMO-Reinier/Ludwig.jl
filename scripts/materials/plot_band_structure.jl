using PlotlyJS
using LaTeXStrings
include("Sr2RuO4.jl")

k_array_x = range(0,0.5,101)
k_array_y = range(0,0.5,101)

band_alpha = zeros(Float64,101,101)

for (i, k_x) in enumerate(k_array_x)
    for (j, k_y) in enumerate(k_array_y)
        band_alpha[i,j] = ham_α([k_x,k_y])
    end
end

#plot(scatter(x=k_array_x, y=band_alpha[:,51], mode="lines"),
#    Layout(xaxis_title=("k_x"), yaxis_title="E(k) (eV)", font=attr(size=18)))

band_gamma_to_X = band_alpha[:,1]
band_X_to_M     = band_alpha[end,:]
band_M_to_gamma = 

total_band = vcat(band_gamma_to_X,band_X_to_M)

plot(scatter(y=total_band, mode="lines"),
    Layout(xaxis_title=("k"), yaxis_title="E(k) (eV)", font=attr(size=18)))