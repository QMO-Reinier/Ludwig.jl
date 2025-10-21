include("electron_electron.jl")
include("/home/lion/Desktop/simulation_code/Code_scaffidi/scripts/materials/Sr2RuO4.jl")
using PlotlyJS
using LaTeXStrings
using LinearAlgebra

# Define parameters
T             = 20  # temperature in Kelvin
n_E_bins      = 4  # number of energy bins
n_theta_bins  = 6 # number of angular bins
annular_width = 6   # 2*annular_width*kbT = width around Fermi surface
material_file = "Sr2RuO4.jl"
data_dir   = "/home/lion/Desktop/simulation_code/Code_scaffidi/data"

bands = [ham_α, ham_β, ham_γ]

outfile = joinpath(@__DIR__, data_dir, "test_data.h5")

main(T, n_E_bins, n_theta_bins, outfile, bands, annular_width)
