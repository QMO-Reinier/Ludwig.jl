using Ludwig
using HDF5
using Interpolations
using StaticArrays
using LinearAlgebra
using ProgressBars
using ProgressMeter

function get_L_imp(T::Real, n_ε::Int, n_θ::Int, outfile::String, bands::Vector{Function}, V_sq, α::Real = 6.0)

    l = Lattices.Lattice([1.0, 0.0], [0.0, 1.0])
    t = @elapsed mesh = Ludwig.bz_mesh(l, bands, kb * T, n_ε, n_θ, 800, α)
    println("Mesh generation took ", t, " seconds.")
    ℓ = length(mesh.patches)
    print(size(mesh.patches))
    print(ℓ)

    # Initialize file - will error if file exists
    h5open(outfile, "cw") do f
        g = create_group(f, "data")
        g["n_ε"] = n_ε
        g["n_θ"] = n_θ
        g["T"] = T
        g["corners"] = copy(transpose(reduce(hcat, mesh.corners)))
        g["momenta"] = copy(transpose(reduce(hcat, map(x -> x.k, mesh.patches))))
        g["velocities"] = copy(transpose(reduce(hcat, map(x -> x.v, mesh.patches))))
        g["energies"] = map(x -> x.e, mesh.patches) 
        g["dVs"] = map(x -> x.dV, mesh.patches)
        #g["corner_ids"] = copy(transpose(reduce(hcat, map(x -> x.corners, mesh.patches))))
        g["corner_ids"] = copy(transpose(reduce(hcat, mesh.corners)))
    end

    T = kb * T # Convert K to eV

    L_imp = zeros(Float64, ℓ, ℓ)
    Ludwig.electron_impurity!(L_imp, mesh.patches, V_sq)

    # Write scattering operator out to file
    h5open(outfile, "cw") do fid
        g = fid["data"]
        g["L_imp"] = L_imp
    end
end
