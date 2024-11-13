# Tracking galaxy evolution with stars Pipeline for TNG50 Simulation

This pipeline processes galaxy merger trees from the TNG50 cosmological simulation, focusing on star particle matching across snapshots to trace properties of central halos. The main HDF5 output files contain galaxy properties over time, enabling detailed analysis of star formation history, mass evolution, and chemical enrichment within the merger tree framework.

### Paths

Please adjust machine-specific paths as necessary, including paths to TNG simulation data and Excel files. Future versions may streamline these paths further, but for now, they are set to work within typical HPC environments.

### Getting Started

The pipeline consists of two main stages, each producing intermediate and final outputs:

1. **Generating Star Particle Matches**: The file `connect_snap.py` identifies galaxies across snapshots by matching star particles to those in a reference galaxy at `snapnum=99`. This produces an Excel file (`MRI_connect_snap_data.xlsx`) that records the galaxy ID with the most star particle matches for each snapshot, creating a reconstructed merger tree.

2. **Extracting and Storing Galaxy Properties**: Based on the reconstructed merger tree, `mergertree_reconstruct.py` extracts properties from each snapshot for the galaxies in the merger tree and stores them in an HDF5 file. This file includes galaxy properties such as mass, metallicity, and star formation rate, tracked over time for each galaxy in the tree.

### Pipeline Steps

The pipeline is organized as follows:

1. **connect_snap.py**  
   - This script matches star particles across snapshots to identify galaxies in each snapshot with the closest relation to the reference galaxy at `snapnum=99`.  
   - Output: `MRI_connect_snap_data.xlsx`, which contains the most particle-matched galaxy IDs for each snapshot.

2. **MRI_mergertree_properties.py**  
   - Reads `MRI_connect_snap_data.xlsx` to retrieve each galaxy’s ID across snapshots.
   - Extracts properties for each galaxy from the TNG50 data, tracking their evolution over time.
   - Output: `MRI_mergertree_evo_properties.hdf5`, an HDF5 file storing various properties for each snapshot, including stellar mass, gas mass, star formation rate, specific star formation rate, and metallicity.

### Running the Pipeline

The entire pipeline can be executed with a job submission system, such as SLURM, by setting up job scripts. Here are the recommended steps:

1. **Extract Star Particle Matches**  
   - Use `connect_snap.py` to create `MRI_connect_snap_data.xlsx`.
   - Run this step in parallel, if possible, by submitting `job-connect.sh` via SLURM:
     ```bash
     sbatch job-connect.sh
     ```

2. **Extract and Store Galaxy Properties**  
   - With `MRI_connect_snap_data.xlsx` in place, run `mergertree_reconstruct.py` to generate the final HDF5 file.
   - Submit this job with `job-properties.sh`:
     ```bash
     sbatch job-properties.sh
     ```

These HDF5 and Excel outputs serve as the foundation for downstream analyses, enabling detailed studies of galaxy evolution within the reconstructed merger trees.
