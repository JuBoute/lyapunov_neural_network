# lyapunov_neural_network
Contain all the relevant code (in python 3) for the paper "Dynamical properties of self-sustained and driven neural networks"
All the necessary libraries are imported at the begining of each files.

## Simulating the networks
The files [Simulation_Adex_driven.py](/Simulation_Adex_driven.py) and [Simulation_Adex_sustained.py](/Simulation_Adex_sustained.py) respectively show the simulation of Driven and Self-sustained Adex networks. 

Both are very similar : they use brian2 to do the simulation, use two different population of neurons : 8000 excitatory, 2000 inhibitory (that have different parameters but identical equations), each neurons have a 5% chance of being conencted to any other, the simulation last for 20s, and we record both the average firing rates and the average adaptation variable of the excitatory neurons (as the inhibitory don't have one). 
The simulations can be done with various seeds that change both the external drive (while it will be statistically identical) and the connectivity (although the same 5% rule is respected). 

The difference lies in the choice of parameters, and in a change in the "external drive" part, as driven networks have a constant one while the self-sustained networks only have a short initial kick.

It is also important to note that not all seeds will lead to self-sustained activity in the self-sustained file. Some suggestions are given, but there would need to repeat simulations for various seeds to fine more (to do so, a simulation time of 20s is not necessary, do not hesitate to change "TotTime")


## Record of the simulation
We recorded two driven simulations with different seeds to show the difference, and one self-sustained. They are respectively in [Simulation_Adex_1](/Simulation_Adex_1) and [Simulation_Adex_2](/Simulation_Adex_2) for the driven networks, and [Simulation_Adex_sustained_1](/Simulation_Adex_sustained_1) for the self-sustained.
As said previously, those files contain the average firing rates and adaptation through time, and can be open with pickle. 

## Analysis of the simulation
The main part of this work was the analysis, and especially how to compute the first lyapunov exponent (FLE) from the average values of the networks. 
The file [Analysis_Lyapunov_Exponent_code.ipynb](/Analysis_Lyapunov_Exponent_code.ipynb) is a Jupyter notebook that summarize this analysis. 

After importing what was needed, we show the functions we will use to compute the FLE.
Then, we show how the analysis is done, first on the 2 driven networks, and then on the self-sustained one. 

This code can be use with other generated files as long as they follow what was done in the Simulation_Adex_.py files.
