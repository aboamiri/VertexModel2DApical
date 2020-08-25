# VertexModel2DApical
A 2D vertex model written in python, for studying vertex stability.
Note: The code is still under development, and may contain bugs that needs to be corrected over time.

runcard.py can be used to run a few simulation types: 
usage: runcard.py sim_type Nxy gamma n_steps dt

for example to run simulation with 8x8 cells, line tension gamma=0.01, for 10000 steps, dt=0.1, run the following:
python3 runcard.py randomized_tension 8 0.01 10000 0.1

Run "python runcard.py -h" for more information on each argument.
