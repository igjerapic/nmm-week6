import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import Counter
from MDAnalysis.core.universe import Universe
from MDAnalysis.lib.distances import calc_angles
from cycler import cycler

plt.style.use('../../scripts/default.mplstyle')
plt.rcParams['axes.prop_cycle'] = plt.cycler(cycler(color = ['#332288', 
                                    '#CC6677',
                                    '#88CCEE',
                                    '#DDCC77', 
                                    '#117733', 
                                    '#882255', 
                                    '#44AA99', 
                                    '#999933', 
                                    '#AA4499',
                                    '#DDDDDD'
                                ]))

def load_universe(data_file, traj_file, format="LAMMPS"):
    return Universe(data_file, traj_file, format=format)


def compute_angles(universe, start=0, end=None):
    angles = []
    residues = universe.residues
    for ts in universe.trajectory[start:end]:
        for i in range(1, len(residues) - 2):
            com1 = residues[i].atoms.center_of_mass(unwrap=True)
            com2 = residues[i + 1].atoms.center_of_mass(unwrap=True)
            com3 = residues[i + 2].atoms.center_of_mass(unwrap=True)
            angle = calc_angles(
                com1[np.newaxis, :],
                com2[np.newaxis, :],
                com3[np.newaxis, :],
                box=universe.dimensions
            )[0]
            angles.append(angle)
    return angles


def process_angle_distribution(angles_rad, bin_step=1):
    angles_deg = [round(np.rad2deg(a), 0) for a in angles_rad]
    counts = Counter(angles_deg)
    x = np.array(sorted(counts.keys()))
    y = np.array([counts[val] / len(angles_deg) for val in x])
    return x, y


def fill_distribution(x, y, step=1, angle_range=(0, 180)):
    x_full = np.arange(angle_range[0], angle_range[1] + step, step)
    y_full = np.zeros_like(x_full, dtype=float)
    indices = np.searchsorted(x_full, x)
    y_full[indices] = y
    return x_full, y_full


def fitting_function(x, a, x0, a1, a2, a3):
    x_rad = (x * np.pi) / 180
    return a * np.exp(a1 * (x_rad - x0) ** 2 + a2 * (x_rad - x0) ** 3 + a3 * (x_rad - x0) ** 4)


def fit_distribution(x, y, initial_guess=[1, 2.0, 1, 1, 1]):
    popt, _ = curve_fit(fitting_function, x, y, p0=initial_guess, maxfev=50000)
    return popt


def plot_distribution(x, y, popt, xlabel, ylabel, save_path=None):
    xx = np.linspace(x.min(), x.max(), 100)
    fig, ax = plt.subplots()
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, max(y) * 1.1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.scatter(x, y, s=100, color='b', alpha=0.9, label="AA")
    ax.plot(xx, fitting_function(xx, *popt), lw=2, ls='--', c='r', label='Fitting')
    ax.legend()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=800)
    plt.show()


def save_distribution(path, x, y):
    np.savetxt(path, np.column_stack((x, y)), fmt="%.1f %.6f")


def main():
    
    # File paths
    data_file = "../PLA_CHARMM.data"
    traj_file = "../unwrapped_traj.dcd"

    # Optional save paths
    save_plot_path = "angle_dist_plot.svg" # e.g., "angle_dist_plot.jpg"
    save_data_path = "angle_dist.txt" # e.g., "angle_distribution.txt"

    # Parameters
    start_frame, end_frame = 0, 2000
    bin_step = 1
    angle_range = (0, 180)

    # Load and compute
    u = load_universe(data_file, traj_file)
    angles_rad = compute_angles(u, start=start_frame, end=end_frame)
    x, y = process_angle_distribution(angles_rad, bin_step=bin_step)
    x_full, y_full = fill_distribution(x, y, step=bin_step, angle_range=angle_range)
    
    popt = fit_distribution(x, y)
    
    kT = 0.59616113  # kT = 0.59616113 Kcal/mol at 300 K
    print('angle_style    quartic')
    print('angle_coeff    1 ',end=' ')
    print(*[f"{val:.5f}" for val in [np.rad2deg(popt[1]), -kT*popt[2], -kT*popt[3], -kT*popt[4]]])

    # Plot if desired
    plot_distribution(
        x, y, popt,
        xlabel=r"$\theta_{A-A-A} \; (\mathrm{^\circ})$",
        ylabel=r"$P_{angle}$",
        save_path=save_plot_path
    )

    # Save if path is specified
    if save_data_path:
        save_distribution(save_data_path, x_full, y_full)


if __name__ == "__main__":
    main()
