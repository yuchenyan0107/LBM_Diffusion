title: Implementation in the APES Framework

In distributed computations it is important to avoid any kind of bottleneck,
as the local computing resources in each computing unit are
usually quite limited in comparison to the total computing power.
For large scale computations that aim to use the full system for their
simulation, it therefore is essential that the complete processing chain is
distributable.
Otherwise any step from the mesh generation over the simulation to the
post-processing might prohibit a successful simulation.
The *TreElM* library builds the basis for a completely
integrated framework that allows the distributed processing of all necessary
simulation steps.
To enable a flexible and convenient user configuration, the Aotus library ([[aotus_module]])
is not only used in *TreElM* header files, but also exposed to the
calling applications.
The layout of the overall framework with the central *TreElM* library is
shown by the schematic figure in the [Motivation](motivation.html).

The solvers read the mesh data, using *TreElM* functions and process the
mesh into a solver specific data structure with attached simulation data.
Results are written in the form of restart files, that follow the elemental
design of the *TreElM* mesh by writing the elemental solution following the
same ordering to disk.
With the help of tracking objects it is possible to write subsets of the mesh
with attached simulation data in arbitrary intervals to disk.
The *TreElM* library not only provides the means to identify the halos to
be exchanged between partitions but also provides various communication
pattern that can be used for the actual exchange.
Due to this encapsulation of the communication, the communication layout can
be easily exchanged.
It is even possible to replace the complete parallel paradigm in the
deployment with only few changes to the code.

Aside from these basic functionalities the library provides several auxiliary
routines, that provide for example logging and debugging facilities.
Finally the post processing tool *Harvester* is used to visualize the
results of the simulation.
Since *Harvester* is a stand alone tool, it can be deployed on specialized
machines, that are more suited for visualization.
Additionally, the separated design ensures the independence of the solvers
from any third party libraries, that are typically required for visualization
outputs.
