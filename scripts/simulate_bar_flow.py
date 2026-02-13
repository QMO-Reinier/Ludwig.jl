# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 12:02:55 2026

@author: regter
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
import os
import h5py

#%%

os.chdir(r"C:\Users\regter\Documents\PhD\Theory\Simulating_boltzmann\Ludwig.jl\data\conductivity")

def load_hdf5_to_numpy(filename):
    """
    Opens an HDF5 file, lists all datasets, and loads them into NumPy arrays.
    
    Returns:
        data_dict: dict[str, np.ndarray] mapping dataset paths to numpy arrays
    """
    data_dict = {}

    def visit_datasets(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(f"Found dataset: {name} | shape={obj.shape} | dtype={obj.dtype}")
            data_dict[name] = obj[()]  # Load data into NumPy array

    with h5py.File(filename, 'r') as f:
        f.visititems(visit_datasets)

    print("\n Loaded all datasets into memory.")
    return data_dict

filename = 'sigma_14K_5n_e_20_n_theta_full_matrix.h5'
data     = load_hdf5_to_numpy(filename)
#%%

q_x      = data["data/q_x_array"]
q_y      = data["data/q_y_array"]
sigma_xx = data["data/sigma_q_xx"]
sigma_xy = data["data/sigma_q_xy"]
sigma_yx = data["data/sigma_q_yx"]
sigma_yy = data["data/sigma_q_yy"]

sigma_tensor_q = np.stack([
    np.stack([sigma_xx, sigma_xy], axis=-1),  # shape (Nx, Ny, 2)
    np.stack([sigma_yx, sigma_yy], axis=-1)   # shape (Nx, Ny, 2)
], axis=-2)  # shape (Nx, Ny, 2, 2)

#%%

Nx = len(q_x)          # points in x (periodic direction)
Ny = len(q_y)          # points in y
delta_q_x = q_x[1] - q_x[0] # assume equal spacing
delta_q_y = q_y[1] - q_y[0] # assume equal spacing
Lx = 2*np.pi/delta_q_x         # system length in x
Ly = 2*np.pi/delta_q_y          # height of strip
W = Ly

dx = Lx / Nx
dy = Ly / Ny

x = np.arange(Nx) * dx
y = np.arange(Ny) * dy

QX, QY = np.meshgrid(q_x, q_y, indexing='ij')

# Boundary y indices
yb_indices = [0, Ny-1]
nb = len(yb_indices)   # number of boundaries = 2

lambda_b = 1e2

# External current (uniform)
j0 = 1.0
j_ext = np.zeros((Nx, 2))
j_ext[:,0] = j0

# FFT in x
j_ext_qx = np.fft.fft(j_ext, axis=0)

# Allocate boundary current in q-space
j_boundary_qx = np.zeros((Nx, nb, 2), dtype=complex)

# Solve boundary equation for each qx
for iq in range(Nx):

    # Build 4x4 matrix (2 components × 2 boundaries)
    M = np.eye(2*nb, dtype=complex)

    for b1 in range(nb):
        for b2 in range(nb):
            y1 = yb_indices[b1]
            y2 = yb_indices[b2]

            # integrate over qy to get kernel in mixed (qx, y)
            for iqy in range(Ny):
                K = sigma_tensor_q[iq, iqy] * np.exp(1j * q_y[iqy] * (y1 - y2))
                M[2*b1:2*b1+2, 2*b2:2*b2+2] += lambda_b * K / Ny

    # RHS
    rhs = np.tile(j_ext_qx[iq], nb)

    # Solve linear system
    sol = np.linalg.solve(M, rhs)

    # Store solution
    for b in range(nb):
        j_boundary_qx[iq, b] = sol[2*b:2*b+2]


# Transform boundary currents back to real x
j_boundary = np.real(np.fft.ifft(j_boundary_qx, axis=0))
#%%

# Initialize full current
j = np.zeros((Nx, Ny, 2), dtype=complex)
j[:,:,0] = j0   # external current

# FFT boundary current in x
j_boundary_qx = np.fft.fft(j_boundary, axis=0)

for iq in range(Nx):
    print(iq)
    for iy in range(Ny):
        for b in range(nb):
            yb = yb_indices[b]

            # integrate over qy
            for iqy in range(Ny):
                K = sigma_tensor_q[iq, iqy] * np.exp(1j * q_y[iqy] * (y[iy] - y[yb]))
                j[iq, iy] -= lambda_b * (K @ j_boundary_qx[iq, b]) / Ny

# inverse FFT in x
j = np.real(np.fft.ifft(j, axis=0))

jx = j[:,:,0]
jy = j[:,:,1]


#%%
import matplotlib.pyplot as plt

extent = [0, Lx, 0, Ly]

# ---- Plot j_x ----
plt.figure()
plt.imshow(jx.T, origin='lower', extent=extent, aspect='auto')
plt.colorbar(label='j_x')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Current Density $j_x(x,y)$')
plt.tight_layout()
plt.show()

# ---- Plot j_y ----
plt.figure()
plt.imshow(jy.T, origin='lower', extent=extent, aspect='auto')
plt.colorbar(label='j_y')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Current Density $j_y(x,y)$')
plt.tight_layout()
plt.show()
