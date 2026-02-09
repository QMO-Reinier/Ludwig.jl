include("electron_electron.jl")
include("materials/Sr2RuO4.jl")
using PlotlyJS
using LaTeXStrings
using LinearAlgebra

# Define parameters
T             = 14  # temperature in Kelvin
n_E           = 4   # n_E-1 is number of energy bins
n_theta       = 10   # n_theta-1 is number of angular bins
annular_width = 6   # 2*annular_width*kbT = width around Fermi surface
material_file = "Sr2RuO4.jl"
data_dir   = raw"\Users\regter\Documents\PhD\Theory\Simulating_boltzmann\Ludwig.jl\data"

bands = [ham_α, ham_β, ham_γ]

outfile = joinpath(@__DIR__, data_dir, "test_data_$(T)K_$(n_E)n_e_$(n_theta)_n_theta.h5")

main(T, n_E, n_theta, outfile, bands, annular_width)
