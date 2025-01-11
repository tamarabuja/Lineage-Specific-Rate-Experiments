# Lineage-Specific-Rate-Experiments

The network visualisation file is a basic network configuration demonstrating the logic of the network. Each node is empty (blue), and can be occupied by a species (green). A seed species is initialised into the network and has a chance of speciating or going extinct, specified with the spec_rate and ext_rate parameters. Each frame generated is a snapshot of the network at each time step in the simulation, representing something like 1 million years during which speciation occurs. 

The lineage specific rate experiment file is slightly different. In this simulation, the relationship between speciation (s) and extinction (e) for any species is defined by the user. In this file it is set to s = e + 0.05, representing a linear relationship between the variables. A number of initial species (seed_species) are generated with random s and e values according to the defined relationship. These species then go onto generate other species (or go exinct) and their s and e properties are passed onto their daughter species. The diversity of each of these resulting lineages is tracked, and then visualised by the curve generator file.
