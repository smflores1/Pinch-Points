# Pinch-Points
This repository contains the simulations supporting predictions in our article on pinch point densities for critical lattice models inside polygons: [arxiv.org/abs/1201.6405](https://arxiv.org/abs/1201.6405).

The simulations section of our article explains the code that you will find in this repository. We have tried to name these files in a self-evident way:
* 'FkQ2' means the Q=2 critical FK (random cluster) model, and 'Perc' means percolation. 
* 'Rec' means rectangle, and 'Hex' means hexagon. 
* 'R2' means a rectangle of aspect ratio 2.

A script, ending with 'Scaled', accompanies each simulation and post-processes the data from the latter. If we are working with a rectangle, then it scales the rectangle down to 2 (length) by 1 (height). If we are working with a hexagon, then it crops the hexagon out of a larger rectangle containing it, deforms it into a regular hexagon, and moves its center to the origin of the x-y plane. In either case, the script also normalizes the data by the value at the center of the rectangle or hexagon.
