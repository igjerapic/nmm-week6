import pickle as pkl
from cycler import cycler

import numpy as np
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
    # determining average polymer size
    rgs = np.loadtxt("equilib_rgs.txt", usecols=2)
    print(np.mean(rgs))
    print(np.std(rgs))

    # plotting linear density 
    df = pkl.load(open("post_process_equilib.pkl", "rb"))
    # "bulk_salt_density" : bulk_salt_density,
    #       "y_vals" : y_vals - 50,
    #       "lin_dense_poly" : lin_dense_poly,
    #       "lin_dense_salt" : lin_dense_salt,
    #       "poly_dense_point": [poly_density_P, salt_density_P],
    #       "supernatent_point":[poly_density_S, salt_density_S]
    #       }
    plt.figure(figsize=(8, 2))
    plt.plot(df["y_vals"], df["lin_dense_poly"][0])
    plt.xlabel("$y~(\sigma)$")
    plt.ylabel(r"$\phi_\text{poly}$ $(m/\sigma^3)$")
    plt.tight_layout()
    plt.savefig("equilib_lindesity.svg")
    plt.show()
if __name__=="__main__":
    main()
