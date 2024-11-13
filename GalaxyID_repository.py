"""
This script processes snapshots from the TNG50-1 simulation to extract and filter galaxy data.
The goal is to create an intermediary file, `snap_all.hdf5`, that stores key properties for each snapshot.
Specifically, `snap_all.hdf5` saves the following data for each snapshot:
    - 'GroupFirstSub': the index of the first subhalo in each halo group.
    - 'GroupNsubs': the number of subhalos in each halo group.
    - 'SubhaloMass': the stellar mass of each subhalo.

After creating `snap_all.hdf5`, the script filters this data to identify galaxies with stellar masses
greater than 1e7. These filtered IDs are saved in a new file, `IDs_above_1e7.hdf5`. This file 
stores, for each snapshot, the IDs of galaxies that meet the mass threshold criterion. 

Usage:
1. Ensure that the TNG simulation data path and output paths are correctly specified.
2. Run the script to generate the `snap_all.hdf5` file and then filter and store IDs in `IDs_above_1e7.hdf5`.
"""

# Import necessary libraries
import sys
import numpy as np
import h5py
from tqdm import tqdm
import illustris_python as il  # TNG simulation package for data handling
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

# Set up system path (update paths as needed)
sys.path.append("/opt/conda/lib/python3.6/site-packages/numpy/")

# Base path for TNG simulation output data (update as needed)
basePath = '/path/to/TNG50-1/output'

# Set numpy print options to show full arrays in print statements
np.set_printoptions(threshold=np.inf)

# Define the range of snapshots to process
snapshots = np.linspace(99, 10, 90).astype(int)

# Process each snapshot to extract data and save it to an HDF5 file
for snap in tqdm(snapshots, desc='Processing snapshots'):
    # Load relevant data from the simulation
    GroupFirstSub = il.groupcat.loadHalos(basePath, snap, fields=['GroupFirstSub'])
    GroupNsubs = il.groupcat.loadHalos(basePath, snap, fields=['GroupNsubs'])
    SubhaloMass = il.groupcat.loadSubhalos(basePath, snap, fields=['SubhaloMassType'])[:, 4] * 1e10 / 0.6774 

    # Save data to an HDF5 file
    output_path = '/path/to/output/snap_all.hdf5'  # update path as needed
    with h5py.File(output_path, 'a') as f:
        # Create a group for each snapshot
        grp = f.create_group(str(snap))
        # Save datasets in the group
        grp.create_dataset('GroupFirstSub', data=GroupFirstSub)
        grp.create_dataset('GroupNsubs', data=GroupNsubs)
        grp.create_dataset('SubhaloMass', data=SubhaloMass)

# Process the HDF5 file to filter subhalo IDs based on mass threshold
input_file_path = '/path/to/output/snap_all.hdf5'  # path to source data
output_file_path = '/path/to/output/IDs_above_1e7.hdf5'  # path to filtered output

# Open the source data file and an output file to store filtered results
with h5py.File(input_file_path, 'r') as f, h5py.File(output_file_path, 'w') as output_file:
    # Loop through each snapshot in the source file
    for snap in tqdm(f.keys(), desc='Processing snaps'):
        # Load data for the current snapshot
        GroupFirstSub = f[snap]['GroupFirstSub'][:]
        GroupNsubs = f[snap]['GroupNsubs'][:]
        SubhaloMass = f[snap]['SubhaloMass'][:]
        
        # Initialize a list to hold IDs of subhalos above the mass threshold
        ids_above_threshold = []
        for i in range(len(GroupFirstSub)):
            for j in range(GroupFirstSub[i], GroupFirstSub[i] + GroupNsubs[i]):
                if j == 0:
                    continue
                # Filter for subhalos with mass greater than 1e7
                if SubhaloMass[j] > 1e7:
                    ids_above_threshold.append(j)
        
        # Save the filtered IDs to the output file
        snap_group = output_file.create_group(f'snap_{snap}')
        snap_group.create_dataset('IDs_above_1e7', data=ids_above_threshold)