import pickle as pkl
from cycler import cycler

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

def main():
    N_SALTS = np.arange(2000, 7001, 1000)
    file_name = "post_process.pkl"
    #plotting linear disityies diagram


    fig, ax = plt.subplots(2, figsize=(4,4))
    for i, N_SALT in enumerate(N_SALTS):
        dir_name = f"salt_{N_SALT}/"
        row = pkl.load(open(dir_name + file_name, "rb"))

        ax[0].plot(row["y_vals"], np.mean(row["lin_dense_poly"], axis = 0))
        ax[1].plot(row["y_vals"], np.mean(row["lin_dense_salt"], axis = 0))

    # plt.legend(ncols=2, columnspacing=0.7, handletextpad=0.2)
    plt.xlabel(r"$\phi_\text{poly}$")
    plt.ylabel(r"$\phi_\text{salt}$")
    plt.tight_layout()
    plt.savefig("ass4_phase-diagram.svg")
    plt.show()

if __name__=="__main__":
    main()