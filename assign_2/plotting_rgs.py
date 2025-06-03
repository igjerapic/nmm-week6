import pickle as pkl
from cycler import cycler

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.style.use('..//scripts/default.mplstyle')
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
fit_colors = ["k", "grey"]
def main():
    #TODO Read in avg RGs for Ns: Good solvent is cutoff 1.12246; Bad is cutoff=2.5
    Ns = np.array([16, 20, 32, 64])
    cutoffs = np.array([1.12246, 2.5])

    file_name = "post_process.pkl"

    rgs = np.zeros((2, len(Ns)))
    rgs_err = np.zeros_like(rgs)


    for i, cutoff in enumerate(cutoffs):
        for j, N in enumerate(Ns):
            dir_name = f"epsilonLJ1.0_cutoff{cutoff}_N{N}/"
            df = pkl.load(open(dir_name + file_name, "rb"))
            
            rgs[i, j] = df["avg_rg"]
            rgs_err[i, j] = df["avg_rg_err"]

    # fitting power law
    power_law = lambda x, a, m: a * x ** m
    slopes = np.zeros_like(cutoffs)
    slopes_err = np.zeros_like(slopes)
    for i in range(len(cutoffs)):
        params, cov = curve_fit(power_law, Ns, rgs[i], sigma = rgs_err[i], absolute_sigma = True)
        slopes[i] = params[-1]
        slopes_err[i] = np.diagonal(cov)[-1] ** 0.5

        fit_vals = power_law(Ns, *params)

        # plotting 
        plt.errorbar(Ns, rgs[i], rgs_err[i], marker=".", capsize=2, linestyle="", label = f"$r_c = {cutoffs[i]}\sigma$")
        plt.plot(Ns, fit_vals, color = fit_colors[i], linestyle=":", label = f"$R_g = bN^\\nu$;\t $\\nu={slopes[i]:.1f} \pm {slopes_err[i]:.1f}$")

    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Radius of Gyration ($\sigma$)")
    plt.xlabel("Number of Beads $N$")
    plt.legend()
    plt.tight_layout()
    plt.savefig("rgs_plot.svg")
    plt.show()

if __name__=="__main__":
    main()