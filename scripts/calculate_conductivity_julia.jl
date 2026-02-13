include(joinpath(@__DIR__, "..","..", "Ludwig.jl/src/properties.jl"))
include(joinpath(@__DIR__, "..","..", "Ludwig.jl/src/constants.jl"))
include("materials/Sr2RuO4.jl")
using HDF5
using Ludwig
using IterativeSolvers
using LinearAlgebra
using ProgressMeter
using Statistics

# Open the datafile
datafile_ee = joinpath(@__DIR__, "..", "data/L_ee", "L_ee_14K_5n_e_20_n_theta.h5")
data_ee = h5open(datafile_ee, "r")

# Extract the data (all have units of eV)
L_ee = read(data_ee["data/L"])             
v    = read(data_ee["data/velocities"])    
E    = read(data_ee["data/energies"])
dV   = read(data_ee["data/dVs"])
T    = read(data_ee["data/T"])

datafile_imp = joinpath(@__DIR__, "..", "data/L_imp", "L_imp_14K_5n_e_20_n_theta.h5")
data_imp = h5open(datafile_imp, "r")

L_imp = read(data_imp["data/L_imp"])

L_tot = L_ee + L_imp
#L_tot = L_ee

#sigma_ee = electrical_conductivity(L_ee, v, E, dV, T*kb, 0.0, [0.0, 0.0])
#sigma_imp = electrical_conductivity(L_imp, v, E, dV, T*kb, 0.0, [0.0, 0.0])
#sigma_tot = electrical_conductivity(L_tot, v, E, dV, T*kb, 0.0, [0.0, 0.0])

# Below we calculate the conductivity as a function of q

N               = 21                     # total number of points
q_x_array_SI    = range(0,100e6,N)     # in m^-1
q_y_array_SI    = range(0,100e6,N)
sigma_q_100_xx  = zeros(ComplexF64,N)    # create empty array to fill up later
sigma_q_110_xx  = zeros(ComplexF64,N)
sigma_q_100_xy  = zeros(ComplexF64,N) 
sigma_q_110_xy  = zeros(ComplexF64,N)
sigma_q_100_yx  = zeros(ComplexF64,N) 
sigma_q_110_yx  = zeros(ComplexF64,N)
sigma_q_100_yy  = zeros(ComplexF64,N) 
sigma_q_110_yy  = zeros(ComplexF64,N)

# Convert to correct units for in electrical_conductivity function
q_x_array   = q_x_array_SI*a/2π
q_y_array   = q_y_array_SI*a/2π

p = Progress(2*N)
for (i,q_x) in enumerate(q_x_array)
    next!(p)
    #sigma_q_100_xx[i] = electrical_conductivity(L_ee, v, E, dV, T*kb, 0.0, [q_x, 0.0])[1,1]
    #igma_q_110_xx[i] = electrical_conductivity(L_ee, v, E, dV, T*kb, 0.0, [q_x/sqrt(2), q_x/sqrt(2)])[1,1]
    #sigma_q_100_xx[i] = electrical_conductivity(L_imp, v, E, dV, T*kb, 0.0, [q_x, 0.0])[1,1]
    #sigma_q_110_xx[i] = electrical_conductivity(L_imp, v, E, dV, T*kb, 0.0, [q_x/sqrt(2), q_x/sqrt(2)])[1,1]
    sigma_q_100_temp  = electrical_conductivity(L_tot, v, E, dV, T*kb, 0.0, [q_x, 0.0])
    next!(p)
    sigma_q_110_temp  = electrical_conductivity(L_tot, v, E, dV, T*kb, 0.0, [q_x/sqrt(2), q_x/sqrt(2)])
    sigma_q_100_xx[i] = sigma_q_100_temp[1,1]
    sigma_q_110_xx[i] = sigma_q_110_temp[1,1]
    sigma_q_100_xy[i] = sigma_q_100_temp[1,2]
    sigma_q_110_xy[i] = sigma_q_110_temp[1,2]
    sigma_q_100_yx[i] = sigma_q_100_temp[2,1]
    sigma_q_110_yx[i] = sigma_q_110_temp[2,1]
    sigma_q_100_yy[i] = sigma_q_100_temp[2,2]
    sigma_q_110_yy[i] = sigma_q_110_temp[2,2]
end

outfile = joinpath(@__DIR__, "..", "data/conductivity", "sigma_14K_5n_e_20_n_theta_no_imp.h5")

h5open(outfile, "cw") do f
        g = create_group(f, "data")
        g["sigma_q_100_xx"] =  sigma_q_100_xx
        g["sigma_q_110_xx"] =  sigma_q_110_xx
        g["sigma_q_100_xy"] =  sigma_q_100_xy
        g["sigma_q_110_xy"] =  sigma_q_110_xy
        g["sigma_q_100_yx"] =  sigma_q_100_yx
        g["sigma_q_110_yx"] =  sigma_q_110_yx
        g["sigma_q_100_yy"] =  sigma_q_100_yy
        g["sigma_q_110_yy"] =  sigma_q_100_yy
        g["q_array"]        =  collect(q_x_array_SI)
end

N               = 21                     # total number of points
q_x_array_SI    = range(0,10π*1e6,N)     # in m^-1, corresponds to 20 um in real space
q_y_array_SI    = range(0,40π*1e6,N)   # corresponds to 5 um in real space
sigma_q_xx      = zeros(ComplexF64, (N,N))
sigma_q_xy      = zeros(ComplexF64, (N,N))
sigma_q_yx      = zeros(ComplexF64, (N,N))
sigma_q_yy      = zeros(ComplexF64, (N,N))

# Convert to correct units for in electrical_conductivity function
q_x_array   = q_x_array_SI*a/2π
q_y_array   = q_y_array_SI*a/2π

p = Progress(N^2)
for (i,q_x) in enumerate(q_x_array)
    for (j,q_y) in enumerate(q_y_array)
        next!(p)
        sigma_q_temp    = electrical_conductivity(L_tot, v, E, dV, T*kb, 0.0, [q_x, q_y])
        sigma_q_xx[i,j] = sigma_q_temp[1,1]
        sigma_q_xy[i,j] = sigma_q_temp[1,2]
        sigma_q_yx[i,j] = sigma_q_temp[2,1]
        sigma_q_yy[i,j] = sigma_q_temp[2,2]
    end
end

outfile = joinpath(@__DIR__, "..", "data/conductivity", "sigma_14K_5n_e_20_n_theta_full_matrix.h5")

h5open(outfile, "cw") do f
        g = create_group(f, "data")
        g["sigma_q_xx"] =  sigma_q_xx
        g["sigma_q_xy"] =  sigma_q_xy
        g["sigma_q_yx"] =  sigma_q_yx
        g["sigma_q_yy"] =  sigma_q_yy
        g["q_x_array"]        =  collect(q_x_array_SI)
        g["q_y_array"]        =  collect(q_y_array_SI)
end
