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
    with open("thermo.pkl", "rb") as f:
        df = pkl.load(f)
    
    keywords = [["Temp"], ['TotEng'], ["E_bond", "E_angle", "E_dihed"]]
    labels = ["Temp", "Press", 'Energy']

    fig = df.plot(x="Step", y=["Temp"], ylabel="Temp", xlabel=r"Step")
    plt.tight_layout()
    plt.savefig('temp_equil.svg')
    plt.show()

    plt.clf()
    fig = df.plot(x="Step", y= ['TotEng'], ylabel="Energy", xlabel=r"Step")
    plt.tight_layout()
    plt.savefig('totEng_equil.svg')
    plt.show()

    plt.clf()
    fig = df.plot(x="Step", y= ["E_bond", "E_angle", "E_dihed"], ylabel="Energy", xlabel=r"Step")
    plt.legend(loc = (0.65, 0.15))
    plt.tight_layout()
    plt.savefig('energies_equil.svg')
    plt.show()

    
    # fig_temp = df.plot(x="Time", y="Temp", ylabel="Temp", figsize=(6,6))
    #plt.savefig('thermo_bondeng.png')
    # plt.show()

    # fig_Press = df.plot(x="Time", y="Press", ylabel="Press")
    # plt.show()

    # fig_energy = df.plot(x="Time", y=['KinEng', 'E_pair'], ylabel="Energy in Reduced Units")
    # plt.show()

if __name__ == '__main__':
    main()