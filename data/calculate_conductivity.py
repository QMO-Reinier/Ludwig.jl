import numpy as np
from matplotlib import pyplot as plt
import h5py
import os

h_bar = 6.50e-16 # Reduced Planck constant in eVs
a     = 3.90e-10 # Lattice constant in m

#%%

os.chdir("P:\Phd\GitHub\Hydro\cuprates_transport\Ludwig.jl\data")

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

filename = 'test_data_14K_5n_e_30_n_theta.h5'
data     = load_hdf5_to_numpy(filename)

#%% In this cell the conductivity is calculated from the collision operator

e_charge = 1.602e-19  # electron charge in Coulomb

L_ee = data["data/L"]
L_ee = L_ee
v = np.vstack((data["data/velocities"][0], data["data/velocities"][1]))

v_SI = v*(1/h_bar)*a # Convert the velocities to units in m/s

sigma_tensor = np.zeros((2,2))
L_ee_in = np.linalg.inv(L_ee)

for i in range(2):
    for j in range(2):
        sigma_tensor[i,j] = 2*e_charge**2 * np.matmul(v_SI[i],np.matmul(L_ee_in,v_SI[j]))
        
print(sigma_tensor)

#%%

def non_local_sigma(L, v, q_vec):
    '''
    This function calculates the non-local conductivity for a given momentum q_vec.
    
    q_vec should be np.array with np.shape(q_vec) = (2,)
    '''
    dim_L = len(L) 
    q_v_mat = np.zeros((dim_L,dim_L))
    
    for i in range(dim_L):
        q_v_mat[i,i] = np.matmul(q_vec,v[:,i])
        
    non_local_sigma_tensor = np.zeros((2,2), dtype=complex)
    operator = np.linalg.inv(L+1j*q_v_mat)
    
    for i in range(2):
        for j in range(2):
            non_local_sigma_tensor[i,j] = 2*e_charge**2 * np.matmul(v[i],np.matmul(operator,v[j]))
    
    return non_local_sigma_tensor

#%% Calculate the non_local conductivity

N = 20
q_max = 130e-6

q_x_array_1 = np.linspace(0,q_max,N)
q_y_array_1 = np.zeros(N)
q_array_1 = np.vstack((q_x_array_1,q_y_array_1))

q_x_array_2 = np.linspace(0,q_max,N)
q_y_array_2 = np.linspace(0,q_max,N)
q_array_2 = np.vstack((q_x_array_2,q_y_array_2))

sigma_xx_array_1 = np.zeros(N, dtype=complex)
sigma_xx_array_2 = np.zeros(N, dtype=complex)

for i in range(len(q_array_1[0])):
    print(i)
    q_vec_1 = q_array_1[:,i]
    q_vec_2 = q_array_2[:,i]
    sigma_xx_array_1[i] = non_local_sigma(L_ee*a/h_bar, v_SI, q_vec_1)[0,0]
    sigma_xx_array_2[i] = non_local_sigma(L_ee*a/h_bar, v_SI, q_vec_2)[0,0]
    

# We normalize the non-local conductivity

norm_sigma_xx_array_1 = sigma_xx_array_1/sigma_xx_array_1[0]
norm_sigma_xx_array_2 = sigma_xx_array_2/sigma_xx_array_2[0]

#%%

plt.figure(figsize=(6,4), dpi=300)
plt.plot(q_x_array_1*1e6, np.real(norm_sigma_xx_array_1), color = "darkgreen", label=r"$\it{q}$ || 100")
plt.plot(q_x_array_2*1e6, np.real(norm_sigma_xx_array_2), color = "mediumblue", label=r"$\it{q}$ || 110")
plt.xlabel(r"$\it{q}$ (a.u.)", fontsize=15)
plt.ylabel(r"$\sigma\it{(q)}/\sigma\it{(0)}$", fontsize=15)
#plt.xticks([0.1,0.6,1.1,], fontsize=15)
#plt.yticks([0.4,0.6,0.8,1], fontsize=15)
plt.legend(frameon=False, fontsize=15)
#plt.xlim((0,1.3))
plt.show()