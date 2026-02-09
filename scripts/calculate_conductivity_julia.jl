include(joinpath(@__DIR__, "..","..", "Ludwig.jl/src/properties.jl"))
include(joinpath(@__DIR__, "..","..", "Ludwig.jl/src/constants.jl"))
include("materials/Sr2RuO4.jl")
using HDF5
using Ludwig
using IterativeSolvers
using LinearAlgebra
using PlotlyJS

# Open the datafile
datafile = joinpath(@__DIR__, "..", "data", "test_data_14K_5n_e_10_n_theta.h5")
data = h5open(datafile, "r")

# Extract the data (all have units of eV)
L_ee = read(data["data/L"])             
v    = read(data["data/velocities"])    
E    = read(data["data/energies"])
dV   = read(data["data/dVs"])
T    = read(data["data/T"])

sigma = electrical_conductivity(L_ee, v, E, dV, T*kb, 0.0, [0.0, 0.0])

# Below we calculate the conductivity as a function of q

N               = 21                   # total number of points
q_x_array_SI    = range(0,10e6,N)       # in m^-1
q_y_array_SI    = range(0,10e6,N)
sigma_q_100_xx  = zeros(ComplexF64,N)    # create empty array to fill up later
sigma_q_110_xx  = zeros(ComplexF64,N)

# Convert to correct units for in electrical_conductivity function
q_x_array   = q_x_array_SI*a/2π
q_y_array   = q_y_array_SI*a/2π


for (i,q_x) in enumerate(q_x_array)
    sigma_q_100_xx[i] = electrical_conductivity(L_ee, v, E, dV, T*kb, 0.0, [q_x, 0.0])[1,1]
end

plt = plot(
    [
        scatter(
            x = q_x_array_SI,
            y = real(sigma_q_100_xx)/real(sigma_q_100_xx[1]),
            mode = "lines",
            name = "σ_xx",
            line = attr(color="green")
        ),
    ],
    Layout(
        xaxis = attr(title="q"),
        yaxis = attr(title="σ_xx"),
        font = attr(size=20),
        #legend=attr(title="Bands"),
    )
)

display(plt)