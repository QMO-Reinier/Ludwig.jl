include("electron_electron.jl")
include("materials/Sr2RuO4.jl")
using PlotlyJS
using LaTeXStrings
using LinearAlgebra

# Define parameters
T             = 20  # temperature in Kelvin
n_E           = 4   # n_E-1 is number of energy bins
n_theta       = 5   # n_theta-1 is number of angular bins
annular_width = 6   # 2*annular_width*kbT = width around Fermi surface
material_file = "Sr2RuO4.jl"
data_dir   = raw"\\file\home_lion$\Regter\Phd\GitHub\Hydro\cuprates_transport\Ludwig.jl\data"

bands = [ham_α, ham_β, ham_γ]

outfile = joinpath(@__DIR__, data_dir, "test_data_$(T)K_$(n_E_bins)n_e_$(n_theta_bins)_n_theta.h5")

main(T, n_E_bins, n_theta_bins, outfile, bands, annular_width)
