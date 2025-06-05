import pickle as pkl

import numpy as np
import beadspring as bsa
from MDAnalysis.analysis import lineardensity
import matplotlib.pyplot as plt


def compute_linear_density(atomgroup, trajectory_slice, axis, binsize = 1.0):
    """
    A fixed version of the compute_linear_density function provided by bsa
    """
    num_frames   = len(trajectory_slice)
    max_box_size = max(atomgroup.universe.dimensions[:3])
    # num of columns is always equal to largest_dimension/binsize
    num_bins = int(max_box_size / binsize)

    density_analyzer = lineardensity.LinearDensity(atomgroup, binsize=binsize)
    density_array = np.zeros((num_frames, num_bins))

    for idx, ts in enumerate(trajectory_slice):
        frame = ts.frame
        density_analyzer.run(frames=[frame], verbose=False)
        density_array[idx,:] = density_analyzer.results[axis]['mass_density']

    return density_array

def main():
    topology = "final.data"
    traj = "traj_salt.lammpstrj"

    u = bsa.setup_universe(topology, traj)
    N_FRAMES = u.trajectory.n_frames

    # Group salts and polymers 
    polymer_atoms = u.select_atoms("type 1 or type 2")
    salt_atoms = u.select_atoms("type 3 or type 4")

    bulk_salt_density = salt_atoms.n_atoms / np.prod(u.dimensions[:3])
    
    # Get linear densities of final N_FRAMES_SLICE frames    
    BINSIZE = 4

    y_vals = np.linspace(min(u.dimensions[:3]), max(u.dimensions[:3]), int(max(u.dimensions[:3])/ BINSIZE))

    lin_dense_poly = compute_linear_density(polymer_atoms, u.trajectory[:1], axis="y", binsize= BINSIZE)
    
    # Defining supernatent (S) and polymer (P) phase
    TOLERANCE = 0.2
    mask_P = lin_dense_poly > (1 - TOLERANCE) * np.max(lin_dense_poly) 
    min_dense_poly = np.min(lin_dense_poly) if np.min(lin_dense_poly) != 0 else TOLERANCE * np.max(lin_dense_poly) 
    mask_S = lin_dense_poly <= min_dense_poly 

    # Extracting polymer density of phases
    poly_density_P = np.mean(lin_dense_poly[mask_P])
    poly_density_S = np.mean(lin_dense_poly[mask_S])

    df = {"bulk_salt_density" : 0.0,
          "y_vals" : y_vals,
          "lin_dense_poly" : lin_dense_poly,
          "lin_dense_salt" : np.zeros_like(lin_dense_poly),
          "poly_dense_point": [poly_density_P, 0],
          "supernatent_point":[poly_density_S, 0]
          }
    
    with open("post_process_equilib.pkl", "wb") as f:
        pkl.dump(df, f)
if __name__ == '__main__':
    main()