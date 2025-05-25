import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.optimize import curve_fit
from MDAnalysis.core.universe import Universe
from MDAnalysis.lib.distances import calc_dihedrals


def load_universe(data_path, traj_path):
    return Universe(data_path, traj_path, format="LAMMPS")


def compute_dihedrals(universe, start=0, end=None):
    angles = []
    for ts in universe.trajectory[start:end]:
        for i in range(1, len(universe.residues) - 3):
            com1 = universe.residues[i - 1].atoms.center_of_mass(unwrap=True)
            com2 = universe.residues[i].atoms.center_of_mass(unwrap=True)
            com3 = universe.residues[i + 1].atoms.center_of_mass(unwrap=True)
            com4 = universe.residues[i + 2].atoms.center_of_mass(unwrap=True)
            angle = calc_dihedrals(
                com1[np.newaxis, :],
                com2[np.newaxis, :],
                com3[np.newaxis, :],
                com4[np.newaxis, :],
                box=universe.dimensions
            )[0]
            angles.append(angle)
    return np.rad2deg(angles)


def build_distribution(data, step=1, range_min=-180, range_max=180):
    rounded = [round(x, 0) for x in data]
    count = Counter(rounded)
    x_vals = np.array(sorted(count.keys()))
    y_vals = np.array([count[x] / len(data) for x in x_vals])

    x_full = np.arange(range_min, range_max + step, step)
    y_full = np.zeros_like(x_full, dtype=float)
    indices = np.searchsorted(x_full, x_vals)
    y_full[indices] = y_vals
    return x_vals, y_vals, x_full, y_full

def fitting_function(x, a0, a1, a2, a3, a4):
    cosx = np.cos(np.radians(x))
    return np.exp(a0 + a1 * cosx + a2 * cosx**2 + a3 * cosx**3 + a4 * cosx**4)

def fitting_function2(x, a0, a1, a2, a3, a4, a5, a6):
    return  np.exp(a0 + a1*np.cos((x*np.pi)/180) + a2*(np.cos((x*np.pi)/180))**2 + a3*(np.cos((x*np.pi)/180))**3 + a4*(np.cos((x*np.pi)/180))**4 + a5*(np.cos((x*np.pi)/180))**5 + a6*(np.cos((x*np.pi)/180))**6) 

def fit_dihedral_function(x, y):
    popt, _ = curve_fit(fitting_function2, x, y, p0=[1,1,10,1,1,1,1], maxfev=50000)
    return popt


def plot_distribution(x, y, model_func, popt, save_path=None):
    xx = np.linspace(x.min(), x.max(), 100)
    fig, ax = plt.subplots(figsize=(5, 3.8))
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, max(y) * 1.1)
    ax.set_xlabel(r"$\Phi_{A-A-A-A} \; (\mathrm{^\circ})$")
    ax.set_ylabel(r"$P_{dihedral}$")
    ax.scatter(x, y, s=100, color='b', alpha=0.9, label="AA")
    ax.plot(xx, model_func(xx, *popt), lw=2, ls='--', c='r', label='Fitting')
    ax.legend()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=800)
    plt.show()


def save_distribution(path, x_full, y_full):
    np.savetxt(path, np.column_stack((x_full, y_full)), fmt="%.1f %.6f")


def main():
    
    data_file = "../PLA_CHARMM.data"
    traj_file = "../unwrapped_traj.dcd"
    
    universe = load_universe(data_file, traj_file)

    start, end = 0, 2000
    dihedrals = compute_dihedrals(universe, start, end)

    x, y, x_full, y_full = build_distribution(dihedrals, step=1)
    
    popt = fit_dihedral_function(x, y)

    # Optional saving
    save_plot_path = "dihedral_dist_plot.svg" # e.g., "dihedral_dist_plot.jpg"
    save_data_path = "dihedral_distribution.txt" # e.g., "dihedral_distribution.txt"

    kT = 0.59616113  # kT = 0.59616113 Kcal/mol at 300 K
    print('dihedral_style    nharmonic')
    print('dihedral_coeff    1 ',end=' ')
    print(*[f"{-kT * p:.5f}" for p in popt])
    
    if save_data_path:
        save_distribution(save_data_path, x_full, y_full)

    plot_distribution(x, y, fitting_function2, popt, save_path=save_plot_path)


if __name__ == "__main__":
    main()
