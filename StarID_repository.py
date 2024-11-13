"""
This script reads from the file `IDs_above_1e7.hdf5`, which contains IDs of galaxies with stellar 
mass greater than 1e7 for each snapshot in a TNG50 simulation. For a specified snapshot, it retrieves 
the ParticleIDs of stars in each qualifying galaxy and saves them in a new HDF5 file. Each galaxy's 
star IDs are stored in a group named after the galaxy ID, under a snapshot-specific file named 
`starID_1e7_gals_snap_<snapshot_number>.hdf5`.

Steps:
1. Set the desired snapshot number (e.g., `snap_20`) and retrieve IDs of galaxies with mass > 1e7.
2. For each galaxy ID, load the star ParticleIDs from the simulation data.
3. Store each galaxy's star ParticleIDs in an output HDF5 file, organized by galaxy ID.

This approach is useful for isolating and storing high-mass galaxy star data from each snapshot,
allowing for efficient access and further analysis.

Requirements:
- Ensure that paths to `IDs_above_1e7.hdf5` and TNG data (`basePath`) are correctly set.
"""

# Import necessary libraries
import h5py
from tqdm import tqdm
import illustris_python as il  # TNG simulation package for data handling

# Define base path for TNG simulation data and input/output files (adjust paths as necessary)
basePath = '/path/to/TNG50-1/output'
input_file_path = '/path/to/IDs_above_1e7.hdf5'
output_file_basepath = '/path/to/StarID/'

# Specify snapshot to process (e.g., 'snap_20')
snapshot_name = 'snap_20'
snapshot_number = int(snapshot_name.split('_')[1])

# Open the HDF5 file containing IDs of galaxies with stellar mass > 1e7
with h5py.File(input_file_path, 'r') as f:
    # Retrieve galaxy IDs for the specified snapshot
    galaxy_ids = f[snapshot_name]['IDs_above_1e7'][:]

    # Initialize progress bar for processing galaxy IDs
    with tqdm(total=len(galaxy_ids), desc=f'Processing snapshot {snapshot_name}') as pbar:
        for subhalo_id in galaxy_ids:
            # Load star ParticleIDs for each galaxy in the specified snapshot
            star_particle_ids = il.snapshot.loadSubhalo(basePath, snapshot_number, subhalo_id, 'star', fields='ParticleIDs')

            # Create an output file for each snapshot and store each galaxy's star ParticleIDs
            output_file_path = f'{output_file_basepath}/starID_1e7_gals_{snapshot_name}.hdf5'
            with h5py.File(output_file_path, 'a') as output_file:
                # Create a group for each galaxy ID within the output file
                galaxy_group = output_file.create_group(str(subhalo_id))
                galaxy_group.create_dataset('StarID', data=star_particle_ids)
            
            # Update the progress bar
            pbar.update(1)