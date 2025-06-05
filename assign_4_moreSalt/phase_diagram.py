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
    N_SALTS = np.arange(2000, 8001, 1000)
    file_name = "post_process.pkl"
    #plotting phase diagram
    salt_conc = np.zeros(2* len(N_SALTS))
    poly_conc = np.zeros_like(salt_conc)
    for i, N_SALT in enumerate(N_SALTS):
        dir_name = f"salt_{N_SALT}/"
        row = pkl.load(open(dir_name + file_name, "rb"))
        poly_dense_point = row["poly_dense_point"]
        supernatent_point = row["supernatent_point"]

        poly_conc[i] = poly_dense_point[0]
        poly_conc[-(i+1)] = supernatent_point[0]

        salt_conc[i] = poly_dense_point[1]
        salt_conc[-(i+1)] = supernatent_point[1]


        plt.scatter(*poly_dense_point, color = f"C{i}", label = f'$\\phi_\\text{{bulk}}={row["bulk_salt_density"]:.2f}$')
        plt.scatter(*supernatent_point, color = f"C{i}")
        plt.plot([poly_dense_point[0], supernatent_point[0]], [poly_dense_point[1],  supernatent_point[1]])
    plt.plot(poly_conc[:len(N_SALTS)], salt_conc[:len(N_SALTS)], ":", color = "grey", zorder=0)
    plt.plot(poly_conc[len(N_SALTS):], salt_conc[len(N_SALTS):], ":", color = "grey", zorder = 0)
    plt.ylim(0.1, 0.65)

    plt.legend(ncols=2, columnspacing=0.7, handletextpad=0.2)
    plt.xlabel(r"$\phi_\text{poly}$ $(m/\sigma^3)$")
    plt.ylabel(r"$\phi_\text{salt}$ $(m/\sigma^3)$")
    plt.tight_layout()
    plt.savefig("ass4_phase-diagram.svg")
    plt.show()

if __name__=="__main__":
    main()