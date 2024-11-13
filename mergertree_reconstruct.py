"""
This script generates an HDF5 file that stores the properties of galaxies across different snapshots 
from a TNG50 simulation. The script is based on an existing Excel file (`MRI_connect_snap_data.xlsx`), 
which stores the galaxy IDs that have the most matching star particles with a selected galaxy at 
`snapnum=99`. The Excel file essentially represents a merger tree that tracks the galaxies 
matching snapnum=99 in each snapshot.

Purpose:
- Using the merger tree data from the Excel file, this script reconstructs a history of 
  each galaxy by calculating and storing its properties (e.g., stellar mass, star formation 
  rate, metallicity) across snapshots in an HDF5 file (`MRI_mergertree_evo_properties.hdf5`).
- The script uses data from different snapshots to extract specific galaxy properties, making 
  it easy to track how these properties evolve over time for each galaxy in the merger tree.

Requirements:
- Paths to `MRI_connect_snap_data.xlsx`, TNG simulation base data, and desired HDF5 output 
  file should be set accordingly.
"""

# Import necessary libraries
import pandas as pd
import h5py
import numpy as np
import gc  # Garbage collection module
from tqdm import tqdm
import illustris_python as il  # TNG simulation package for data handling

# %matplotlib inline
import sys
sys.path.append("/opt/conda/lib/python3.6/site-packages/numpy/")
import numpy as np
import h5py
from tqdm import tqdm
from Project_galaxy import deIllustrisTNG_galaxy
# %matplotlib widget
run = 'TNG50'
basePath = '/home/tnguser/sims.TNG/TNG50-1/output'
# 设置打印选项，将省略号取消
np.set_printoptions(threshold=np.inf)

# Define constants and paths
h = 0.6774
ZH_sun = (0.0196 / 0.7381)
full_snapshots = [99, 91, 84, 78, 72, 67, 59, 50, 40]  # List of snapshots for detailed data extraction

# Replace these with the paths for your environment
excel_file_path = '/path/to/MRI_connect_snap_data.xlsx'  # Path to Excel file with merger tree data
hdf5_file_path = '/path/to/MRI_mergertree_evo_properties.hdf5'  # Path to the output HDF5 file

# Load the Excel file with merger tree data
df = pd.read_excel(excel_file_path)

# Create the HDF5 file to store galaxy properties across snapshots
with h5py.File(hdf5_file_path, 'w') as hdf_file:
    # Iterate over each row in the DataFrame, representing each tracked galaxy
    for index, row in tqdm(df.iterrows(), total=df.shape[0], desc='Processing rows'):
        # Use 'Snap99_ID' as the group name in HDF5 for each galaxy
        group_name = str(int(row['Snap99_ID']))
        group = hdf_file.create_group(group_name)
        
        # Initialize lists for storing galaxy properties across snapshots
        snaps, subhalo_ids, numbers = [], [], []
        Star_mass_all, DM_mass_all, Gas_mass_all, SFR_all, sSFR_all, Re_all, H_all, Z_all, ZH_all = [], [], [], [], [], [], [], [], [], []
        Effective_star_mass_all, sSFR_inner_list, sSFR_outer_list = [], [], []
        Z_inner_gas_list, Z_outer_gas_list, Z_inner_star_list, Z_outer_star_list = [], [], [], []
        Z_Tagged_inner_gas_list, Z_Tagged_outer_gas_list, Z_Tagged_inner_star_list, Z_Tagged_outer_star_list = [], [], [], []
        M200_all, Rvir_all = [], []
        x, y, z, vx, vy, vz = [], [], [], [], [], []

        # Loop through snapshots, tracking galaxy properties
        for snap in range(99, 29, -1):
            snap_column = f'Snap{snap}_ID'
            number_column = f'Number({snap})'

            # Check if the current snapshot is part of the merger tree and has matching data
            if snap_column in df.columns and pd.notna(row[snap_column]):
                # Record snapshot properties for galaxies that match the criteria
                snaps.append(int(snap))
                subhalo_ids.append(int(row[snap_column]))
                numbers.append(int(row[number_column]))

                # Load subhalo and halo properties
                Subhalo = il.groupcat.loadSingle(basePath, int(snap), subhaloID=int(row[snap_column]))
                haloID = Subhalo['SubhaloGrNr']
                halo = il.groupcat.loadSingle(basePath, int(snap), haloID=haloID)
                
                # Calculate various galaxy properties and append to lists
                Star_mass = Subhalo['SubhaloMassInRadType'][4] * 1e10 / h
                Effective_star_mass = Subhalo['SubhaloMassInHalfRadType'][4] * 1e10 / h
                DM_mass = Subhalo['SubhaloMassInRadType'][1] * 1e10 / h
                Gas_mass = Subhalo['SubhaloMassInRadType'][0] * 1e10 / h
                SFR = Subhalo['SubhaloSFR']
                sSFR = np.log10((SFR / Star_mass) / 1e-8)
                Re = Subhalo['SubhaloHalfmassRadType'][4] * 1000 / h
                H = Subhalo['SubhaloStarMetalFractions'][0]
                Z = np.sum(Subhalo['SubhaloStarMetalFractions'][-8:])
                ZH = np.log10((Z / H) / ZH_sun)
                
                # Additional data for selected snapshots in full_snapshots
                if snap in full_snapshots:
                    star = il.snapshot.loadSubhalo(basePath, snap, int(row[snap_column]), 'star', fields=['Coordinates', 'Masses', 'GFM_Metals', 'GFM_Metallicity', 'GFM_MetalsTagged'])
                    gas = il.snapshot.loadSubhalo(basePath, snap, int(row[snap_column]), 'gas', fields=['Coordinates', 'Masses', 'StarFormationRate', 'GFM_Metals', 'GFM_Metallicity', 'GFM_MetalsTagged'])

                    # Additional processing for tagged metals, gas, and star particles
                    # [Code for tagged metal processing and filtering omitted for brevity]
                
                # Record halo-related properties
                Rvir = halo['Group_R_Crit200'] / h
                M200 = halo['Group_M_Crit200'] * 1e10 / h

                # Position and velocity relative to halo center
                Pos = Subhalo['SubhaloPos']
                Vel = Subhalo['SubhaloVel']
                cen_Pos = halo['GroupPos']
                cen_Vel = halo['GroupVel']

                # Calculate distances and velocities relative to halo center
                dist = np.sqrt(np.sum((Pos - cen_Pos)**2))
                dist_normalized = (dist / h) / Rvir

                x.append((Pos[0] - cen_Pos[0]) / h)
                y.append((Pos[1] - cen_Pos[1]) / h)
                z.append((Pos[2] - cen_Pos[2]) / h)
                vx.append(Vel[0] - cen_Vel[0])
                vy.append(Vel[1] - cen_Vel[1])
                vz.append(Vel[2] - cen_Vel[2])

                # Append calculated properties to lists
                Star_mass_all.append(Star_mass)
                Effective_star_mass_all.append(Effective_star_mass)
                DM_mass_all.append(DM_mass)
                Gas_mass_all.append(Gas_mass)
                SFR_all.append(SFR)
                sSFR_all.append(sSFR)
                Re_all.append(Re)
                H_all.append(H)
                Z_all.append(Z)
                ZH_all.append(ZH)
                
                # Additional lists for specific snapshot data processing omitted for brevity

                M200_all.append(M200)
                Rvir_all.append(Rvir)

        # Create datasets for each calculated property within the HDF5 group
        group.create_dataset('Snap', data=np.array(snaps))
        group.create_dataset('SubhaloID', data=np.array(subhalo_ids))
        group.create_dataset('Number', data=np.array(numbers))
        group.create_dataset('Effective_star_mass', data=np.array(Effective_star_mass_all))
        
        # [Additional dataset creations omitted for brevity]

        # Cleanup after each snapshot to manage memory
        del snaps, subhalo_ids, numbers, Star_mass_all, DM_mass_all, Gas_mass_all, SFR_all, sSFR_all, Re_all, H_all, Z_all, ZH_all, x, y, z, vx, vy, vz
        del Effective_star_mass_all, sSFR_inner_list, sSFR_outer_list, Z_inner_gas_list, Z_outer_gas_list, Z_inner_star_list, Z_outer_star_list
        del Z_Tagged_inner_gas_list, Z_Tagged_outer_gas_list, Z_Tagged_inner_star_list, Z_Tagged_outer_star_list, M200_all, Rvir_all
        gc.collect()  # Trigger garbage collection to free memory

print(f'HDF5 file created at {hdf5_file_path}')
