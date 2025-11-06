using PlotlyJS
using LaTeXStrings
using LinearAlgebra
using Interpolations
include("Sr2RuO4.jl")

band = ham_γ

k_array_x = range(0,0.5,101)
k_array_y = range(0,0.5,101)

vel = ForwardDiff.gradient(x -> band(x[1], x[2]), [k_array_x, k_array_y])

#=

N = 1001
x = LinRange(-0.5, 0.5, N)
E = Array{Float64}(undef, N, N)

for i in 1:N, j in 1:N
    E[i, j] = band([x[i], x[j]]) # Get eigenvalues (bands) of each k-point
end

itp = interpolate(E, BSpline(Cubic(Line(OnGrid()))))

=#