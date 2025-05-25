import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from MDAnalysis.core.universe import Universe
from collections import Counter
from MDAnalysis.lib.distances import calc_bonds


def load_universe(data_file, traj_file, format="LAMMPS"):
    return Universe(data_file, traj_file, format=format)


def compute_bond_distances(u, start=0, end=None):
    bond_lengths = []
    residues = u.residues
    for ts in u.trajectory[start:end]:
        for i in range(1, len(residues) - 1):
            res1 = residues[i]
            res2 = residues[i + 1]
            com1 = res1.atoms.center_of_mass(unwrap=True)
            com2 = res2.atoms.center_of_mass(unwrap=True)
            bond = calc_bonds(com1[np.newaxis, :], com2[np.newaxis, :], box=u.dimensions)[0]
            bond_lengths.append(bond)
    return bond_lengths


def process_bond_distribution(bonds, decimal=1):
    rounded = [round(b, decimal) for b in bonds]
    cnt = Counter(rounded)
    x = np.array(sorted(cnt.keys()))
    y = np.array([cnt[val] / len(bonds) for val in x])
    return x, y


def fill_distribution(x, y, x_range, step):
    x_full = np.arange(*x_range, step)
    y_full = np.zeros_like(x_full, dtype=float)
    indices = np.searchsorted(x_full, x)
    y_full[indices] = y
    return x_full, y_full


def fitting_function(x, a, x0, a1, a2, a3):
    class2_bond = a1 * (x - x0) ** 2 + a2 * (x - x0) ** 3 + a3 * (x - x0) ** 4
    return a * np.exp(class2_bond)


def fit_distribution(x, y, initial_guess=[1, 3.3, 1, 1, 1]):
    popt, _ = curve_fit(fitting_function, x, y, p0=initial_guess, maxfev=50000)
    return popt


def plot_distribution(x, y, popt, xlabel, ylabel, output_file=None):
    xx = np.linspace(x.min(), x.max(), 100)
    fig, ax = plt.subplots(figsize=(5, 3.8))
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, max(y) * 1.1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.scatter(x, y, s=100, color='b', alpha=0.9, label="AA")
    ax.plot(xx, fitting_function(xx, *popt), lw=2, ls='--', c='r', label='Fitting')
    ax.legend()
    if output_file:
        plt.savefig(output_file, bbox_inches='tight', dpi=800)
    plt.show()


def save_distribution(file_path, x_full, y_full):
    np.savetxt(file_path, np.column_stack((x_full, y_full)), fmt="%.1f %.6f")


def main():
    
    # Input files
    data_file = "../PLA_CHARMM.data"
    traj_file = "../unwrapped_traj.dcd"

    # Parameters
    start, end = 0, 2000
    step = 0.1
    x_range = (2.5, 5.0)
    save_plot_path = "bond_dist_plot.svg"  # e.g., "bond_dist_plot.jpg"
    save_data_path = "bond_distribution.txt"  # e.g., "bond_distribution.txt"

    u = load_universe(data_file, traj_file)
    bonds = compute_bond_distances(u, start, end)
    x, y = process_bond_distribution(bonds)
    x_full, y_full = fill_distribution(x, y, x_range, step)
    
    popt = fit_distribution(x, y)
    
    plot_distribution(x, y, popt, xlabel=r"$l_{A-A} \; (\mathrm{\AA})$", ylabel=r"$P_{bond}$", output_file=save_plot_path)
    
    kT = 0.59616113  # kT = 0.59616113 Kcal/mol at 300 K
    print('bond_style    class2')
    print('bond_coeff    1 ',end=' ')
    print( *[f"{val:.5f}" for val in [popt[1], -kT*popt[2], -kT*popt[3], -kT*popt[4]]])
    
    if save_data_path:
        save_distribution(save_data_path, x_full, y_full)


if __name__ == "__main__":
    main()
