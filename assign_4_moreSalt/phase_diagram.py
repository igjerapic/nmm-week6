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
    N_SALTS = np.arange(0, 1501, 150)
    file_name = "post_process.pkl"
    df = pd.DataFrame(columns = ["bulk_salt_density" ,
                                 "y_vals",
                                 "lin_dense_poly",
                                 "lin_dense_salt",
                                 "poly_dense_point",
                                 "supernatent_point"] )
    
    df = {}
    #plotting phase diagram
    for i, N_SALT in enumerate(N_SALTS):
        dir_name = f"salt_{N_SALT}/"
        row = pkl.load(open(dir_name + file_name, "rb"))
        poly_dense_point = row["poly_dense_point"]
        supernatent_point = row["supernatent_point"]
        plt.scatter(*poly_dense_point, color = f"C{i}")
        plt.scatter(*supernatent_point, color = f"C{i}")

    plt.show()

if __name__=="__main__":
    main()