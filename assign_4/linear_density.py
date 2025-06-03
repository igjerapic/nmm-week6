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
        density_analyzer.run(frames=[idx], verbose=False)
        density_array[idx,:] = density_analyzer.results[axis]['mass_density']

    return density_array

def main():
    # Define the topology and trajectory files

    topology = "final.data"
    traj = "traj_salt.lammpstrj"

    u = bsa.setup_universe(topology, traj)


    # Group salts and polymers 
    polymer_atoms = u.select_atoms("type 1 or type 2")
    salt_atoms = u.select_atoms("type 3 or type 4")

    # Get linear densities of final N_FRAMES_SLICE frames    
    N_FRAMES_SLICE = 1
    lin_dense_poly = compute_linear_density(polymer_atoms, u.trajectory[-N_FRAMES_SLICE:], axis="y", binsize= 1.0)
    lin_dense_salt = compute_linear_density(salt_atoms, u.trajectory[-N_FRAMES_SLICE:], axis="y", binsize= 1.0)

    # Defining supernatent (S) and polymer (P) phase
    tol = 0.2
    mask_P = lin_dense_poly > (1 - tol) * np.max(lin_dense_poly) 
    min_dense_poly = np.min(lin_dense_poly) if np.min(lin_dense_poly) != 0 else tol * np.max(lin_dense_poly) 
    mask_S = lin_dense_poly <= min_dense_poly 
    print(min_dense_poly)
    print(np.partition(lin_dense_poly, 2))

    # Extracting salt and polymer density of phases
    poly_density_P = np.mean(lin_dense_poly[mask_P])
    salt_density_P = np.mean(lin_dense_salt[mask_P])

    poly_density_S = np.mean(lin_dense_poly[mask_S])
    salt_density_S = np.mean(lin_dense_salt[mask_S])

    print("Polymer rich:", [poly_density_P, salt_density_P])
    print("Supernatent:", [poly_density_S, salt_density_S])


    # # Saving to pickle file
    # df_post_process = { "time": time,
    #                     "rgs": rgs,
    #                     "avg_rg": avg_rg,
    #                     "avg_rg_err": avg_rg_err,
    #                     }
    
    # with open("post_process.pkl", "wb") as f:
    #     pkl.dump(df_post_process, f)

if __name__ == '__main__':
    main()