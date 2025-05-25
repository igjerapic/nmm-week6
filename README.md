# nmm-week6

Colloids are fun, but let's go up in scale and model complexity! This assignment is similar in scope to the one of last week, but we focus on polymers with coarse grained models.

## Assignment 1 - Who needs Atoms?

Develop a coarse grained model for a simple polymer (polylactic acid) from a higher-resolution all-atom model. The complex macromolecular structure means that the chemical connectivity must be considered by using multiple beads with bonded interactions. For both bonded and non-bonded interactions, we follow a coarse-graining scheme called (Iterative) Boltzmann Inversion. With number of iterations N=1... Credit to Yang Wang (wangyang.cgmd@gmail.com) for the original scripts.

### Instructions

1a. Load the `PLA_CHARMM.data` atomistic structure in Ovito and inspect it. How many atoms belong to the same monomer? How many monomers are present in the entire chain?

1b. Define the bonded interactions of the CG model.

(i) Run the script PLA_CHARMM.in of a single atomistic chain.

(ii) From the sampled configurations, use the scripts in the `prob_distr_calc/` directory to plot the distributions of all bonds, angles, and dihedrals between the centers of mass of the CG beads. The scripts plot the distributions for the "effective" CG beads from the all-atom simulations, fit them to an analytical function and output the fitting parameters. How does MDanalysis identify the CG beads? Tip: look into the compute_bond_distances() function in the compute_bonds.py file.

(iii) From your fits, write the equations used for the bonded potential terms of the CG force field, specifying the values of all parameters obtained. Discuss briefly also what you would have to do to then calibrate the non-bonded interactions.

1d (OPTIONAL, EASY BUT TAKES SOME TIME). Calibrate the non-bonded interactions to complete the CG force field. Run both AA and CG simulations of several chains (density ~0.9, temperature 300K), and compare the distributions of bonds, angles, and radial distribution function between the AA and the CG models.

## Assignment 2 - Walking Randomly

In this assignment we cover the very basics of polymer physics, focusing on the conformation in solution of polymers with varying architecture and molecular interactions. For this exercise, we introduce the bead-spring model for polymers.

### Instructions

2a. Start from the `run_single_chain.in` file and complete the script with the values of all relevant missing parameters. Obtain the scaling law for the [radius of gyration](https://library.fiveable.me/physical-chemistry-ii/unit-7/polymer-conformations-radius-gyration/study-guide/ZZkPanXFxDIZiVtt) of a bead-spring polymer chain in solution with varying length N. Report two plots, one for good solvent and one for poor solvent conditions. Compare your results to theoretical models.

To build the simulation box with one chain, first create a linear chain data and molecule files using `create_polymer.py` script. Set the sidechain length (SC_length) = 0. The data file is useful to inspect the created chain in Ovito. Execute `run_lammps.sh` to run the `pack_chains.in` script that will pack the molecule into a larger box. Make sure the name of the molecule file is correct in the `pack_chains.in` script.

## Assignment 3 - Digital Breaking Bad

Thanks to advanced synthesis techniques, it is possible to design polymers with tunable architecture and chemical structure, so-called block copolymers, that can assemble into highly-controlled structures. Our synthesis technique is the Python script used in assignment 2, let's invent new polymers that chemists can only dream of and let's see what they assemble into!

### Instructions

3a. Explore in more detail the script `create_polymer.py`.

(i) What are the tunable parameters of the script? Which polymer architectures can be generated?

(ii) Using the `create_polymer.py` script and `pack_chains.in`, generate at least three different data files, each containing a mixture of block copolymers of your choice, with tunable non-bonded interactions for each bead type. Make snapshots with Ovito and discuss the self-assembly behavior you expect.

3b. Choose two of your structures. Adapt the `run.in` file to define all needed non-bonded parameters of the force field and extend the run time. Run the simulation to study their self-assembly behavior in solution. Report your results both with snapshots and quantitative metrics. Include proof that your system is equilibrated.

Pro tips: make sure to have enough polymers to get good statistics, but not so many that the equilibration takes forever. Use the literature to guide your choice. Stay away from high-density phases if you don't want to wait forever to reach an equilibrium phase. Adjust the box length and number of chains in the `pack_chains.in` script to get the desired system. Before running long full simulations, test a much shorter version of your script to make sure it works as intended.


## Assignment 4 - Opposites Attract

Assignment 3 dealt with so-called "simple coacervation", the assembly of same-type macromolecules/monomers. Complex coacervation, where two (or more) different macromolecules assemble in water is also ubiquitous in materials and biological phenomena. [It is typically observed between polymers with oppositely charged groups, and has received attention lately both from fundamental and application perspectives](https://doi.org/10.1039/D0SM00001A). This can be observed with simple models of polyelectrolytes in water, which you will do in this assignment that uses a slab setup.

### Instructions

4a. As you saw in assignment 3, it can take forever to achieve complete phase separation. In the case of complex coacervation, it is even slower due to the Coulomb calculations (recap question: how does the computational cost scale in the presence of charges?). So here, we start from a system that is already phase separated and measure the critical concentration of salt ions necessary to dissolve the polymer-rich phase. See the lecture slides for a schematic of what the binodal curve looks like.

(i) Inspect the `equilibrated.data` file and determine the number of polymers of each type, the number of monomers in one chain, and the magnitude of the charge on the two types of polymers.

(ii) Measure the relative densities of the polymer-rich phase and the solvent-rich phase. You can either use the linear_density function in the bead-spring module or use the histogram modifier in Ovito to calculate the densities. Measure the radius of gyration of the polymers and discuss potential finite-size effects of your simulation setup.

4b. Time to dissolve the coacervate. Repeat the simulations by adding an increasing amount of salt, simulated in the form of individual beads with charge +1 and -1 (always in equal amounts, to ensure charge neutrality of the system). Use the `add_salt.in` script to add salt to `equilibrated.data`. For each simulation, report in the same plot the densities of the polymer-rich phase and the solvent-rich phase. Doing so will allow you to build the binodal curve for the complex coacervation process. From this plot, estimate the critical salt concentration at which the complex coacervation does not happen anymore. Compare the salt concentration found to experimentally reported values and discuss potential limitations of this model.
