import pickle as pkl

import numpy as np
import beadspring as bsa
import matplotlib.pyplot as plt

def main():
    # Define the topology and trajectory files

    topology = "final.data"
    traj_lin = "traj_lin.lammpstrj"

    u = bsa.setup_universe(topology, traj_lin)

    N_FRAMES = u.trajectory.n_frames
    N_ATOMS = u.atoms.n_atoms

    # Initialise the position and time arrays
    positions = np.zeros((N_FRAMES, N_ATOMS, 3))
    time = np.zeros(N_FRAMES)
    rgs = np.zeros(N_FRAMES)
    # Loop over the trajectory and load the positions
    for i,traj in enumerate(u.trajectory):                          
        positions[i] = u.atoms.positions   
        time[i] = u.trajectory.ts.data['time']

        # calculated Radisu of gyration from eigvals of  gyration tensor
        _, eigvals =  bsa.compute_gyration_tensor(positions[i])
        rgs[i] = bsa.calculate_rg2(*eigvals) ** 0.5

    # averaging over trajectory 
    avg_rg = np.mean(rgs)
    avg_rg_err = np.std(rgs)

    # Saving to pickle file
    df_post_process = { "time": time,
                        "rgs": rgs,
                        "avg_rg": avg_rg,
                        "avg_rg_err": avg_rg_err,
                        }
    
    with open("post_process.pkl", "wb") as f:
        pkl.dump(df_post_process, f)

if __name__ == '__main__':
    main()
