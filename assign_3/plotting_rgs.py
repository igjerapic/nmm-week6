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
    file_name = "rgs.txt"
    dirs = ["grafted", "linear"]

    rgs = {}
    rgs_err= {}
    for dir in dirs:
        rgs_vals = np.loadtxt(f"{dir}/" + file_name, usecols=(2)).T
        rgs[dir] = np.mean(rgs_vals)
        rgs_err[dir] = np.std(rgs_vals)
    print(rgs)
    print(rgs_err)

if __name__=="__main__":
    main()