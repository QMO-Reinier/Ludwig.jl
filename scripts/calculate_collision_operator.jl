include("electron_electron.jl")
include("electron_impurity.jl")
include("materials/Sr2RuO4.jl")
using PlotlyJS
using LaTeXStrings
using LinearAlgebra

# Define parameters
T             = 14  # temperature in Kelvin
n_E           = 4   # n_E-1 is number of energy bins
n_theta       = 5   # n_theta-1 is number of angular bins
annular_width = 6   # 2*annular_width*kbT = width around Fermi surface
material_file = "Sr2RuO4.jl"
data_dir_L_ee   = raw"\Users\regter\Documents\PhD\Theory\Simulating_boltzmann\Ludwig.jl\data\L_ee"
data_dir_L_imp   = raw"\Users\regter\Documents\PhD\Theory\Simulating_boltzmann\Ludwig.jl\data\L_imp"

bands = [ham_α, ham_β, ham_γ]

outfile_L_ee = joinpath(@__DIR__, data_dir_L_ee, "L_ee_$(T)K_$(n_E)n_e_$(n_theta)_n_theta.h5")
outfile_L_imp = joinpath(@__DIR__, data_dir_L_imp, "L_imp_$(T)K_$(n_E)n_e_$(n_theta)_n_theta.h5")

main(T, n_E, n_theta, outfile_L_ee, bands, annular_width)

function V_squared(a,b)
    return ...
end

get_L_imp(T, n_E, n_theta, outfile_L_imp, bands, V_squared, annular_width)

