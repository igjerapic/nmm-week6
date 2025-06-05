import pickle as pkl
from cycler import cycler


import matplotlib.pyplot as plt
plt.style.use('../scripts/default.mplstyle')

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
    file = "thermo.pkl"
    dirs = ["grafted", "linear"]
    with open("grafted/thermo.pkl", "rb") as f:
        df_bottle = pkl.load(f)
    with open("linear/thermo.pkl", "rb") as f:
        df_linear = pkl.load(f)

    print(df_linear.keys())
    keywords = [["Temp"], ["Press"], ['TotEng'], ["Density"]]
    labels = ["Temp ($T k_B/\epsilon$)", "Press ($P\sigma^3/\epsilon$)", 'Energy ($E/\epsilon$)', "Density ($\rho \sigma^3/ m$)"]

    for y, ylabel in zip(keywords, labels):
        plt.plot(df_linear["Step"], df_linear[y], label = "linear")
        plt.plot(df_bottle["Step"], df_bottle[y], label = "grafted")
        plt.ylabel(ylabel)
        plt.xlabel("Step")
        plt.tight_layout()
        plt.legend()
        plt.savefig(f"{y[0]}_equilib.svg")
        plt.clf()
        #plt.show()
    # fig_temp = df.plot(x="Time", y="Temp", ylabel="Temp", figsize=(6,6))
    #plt.savefig('thermo_bondeng.png')
    # plt.show()

    # fig_Press = df.plot(x="Time", y="Press", ylabel="Press")
    # plt.show()

    # fig_energy = df.plot(x="Time", y=['KinEng', 'E_pair'], ylabel="Energy in Reduced Units")
    # plt.show()

if __name__ == '__main__':
    main()