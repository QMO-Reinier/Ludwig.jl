import numpy as np
from matplotlib import pyplot as plt
import h5py
import os

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

    print("\n✅ Loaded all datasets into memory.")
    return data_dict

filename = 'test_data_20K_4n_e_5_n_theta.h5'
data     = load_hdf5_to_numpy(filename)

#%%

e_charge = 1.602e-19  # electron charge in Coulomb

L_ee = data["data/L"]
v = np.vstack((data["data/velocities"][0], data["data/velocities"][1]))

sigma_tensor = np.zeros((2,2))

for i in range(2):
    for j in range(2):
        L_ee_in = np.linalg.inv(L_ee)
        sigma_tensor[i,j] = e_charge**2 * np.matmul(v[i],np.matmul(L_ee_in,v[j]))
        
print(sigma_tensor)

#%%

def non_local_sigma(L, v_vec, q_vec, v):
    '''
    This function calculates the non-local conductivity as a function of 
    momentum (q).
    v_vec and q_vec must be same dimension
    '''
    non_local_sigma_tensor = np.zeros((2,2,len(v_vec)))
    
    for i in range(2):
        for j in range(2):
            operator = np.linalg.inv(L_ee+1j*np.matmul(q_vec,v_vec))
            non_local_sigma_tensor[i,j] = e_charge**2 * np.matmul(v[i],np.matmul(L_ee_in,v[j]))
    
    return non_local_sigma_tensor

#%%

