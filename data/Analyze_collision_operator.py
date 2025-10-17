import numpy as np
from matplotlib import pyplot as plt
import h5py
import os

os.chdir("/home/lion/Desktop/simulation_code/Code_scaffidi/data")

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

filename = 'Collision_operator_20K_Sr2RuO4_first.h5'
data     = load_hdf5_to_numpy(filename)

#%%

L_ee = data["data/L"]
v_x  = data["data/velocities"][0]
v_y  = data["data/velocities"][0]

sigma_tensor = np.zeros((2,2))

for i in range(2):
    for j in range(2):
        sigma_tensor[i,j] = 
