include(joinpath(@__DIR__, "..","..", "Ludwig.jl/src/properties.jl"))
include(joinpath(@__DIR__, "..","..", "Ludwig.jl/src/constants.jl"))
using HDF5
using Ludwig

# Open the datafile
datafile = joinpath(@__DIR__, "..", "data", "test_data_14K_5n_e_20_n_theta.h5")
data = h5open(datafile, "r")

# Extract the data
L_ee = read(data["data/L"])
v    = read(data["data/velocities"])
E    = read(data["data/energies"])
dV   = read(data["data/dVs"])
T    = read(data["data/T"])

# Convert to SI units
v = 

electrical_conductivity(L, v, E, dV, T, ω = 0.0, q = [0.0, 0.0])