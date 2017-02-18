# Moleidoscope (Molecular Kaleidoscope)

## About

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://github.com/kbsezginel/Moleidoscope/blob/master/doc/Figures/Fig3.gif">

A key challenge in chemistry is to design molecules with given shapes and sizes comprising many components. One approach to designing structures with such a large number of components is to leverage the typically highly symmetrical nature of supramolecular structures (Figure 1). “Moleidoscope” interactively applies point group symmetry operations to generate hypothetical supramolecular structures in silico. Starting with a simple organic compound, the molecule is replicated using mirror operations, and by selecting multiple mirrors (3D planes), prisms can be formed to have multiple copies of the molecule oriented symmetrically. Moreover, as is the case with kaleidoscope, by rotating these mirrors even more copies can be obtained with many different symmetries. 
***
![alt text][Fig1]

**Figure 1:** Supramolecular structure generation process. First a library of supramolecular structures are decomposed into their building blocks which are simplified as building blocks. Using this library of building blocks new structures can be discovered by assembling these building blocks in different ways.
***
“Moleidoscope” performs these symmetry operations by selecting different organic compounds that can be found in [HostDesigner][HD], a software developed by Dr. Benjamin Hay in Oak Ridge National Laboratory. Also, it acts as a Python wrapper for HostDesigner which can be used to find linking fragments between the arrangement of molecules. In this way, different hypothetical supramolecular structures can be discovered.

### Supramolecular Cages

![alt text][Fig2]

Using _moleidoscope_ supramolecular cages with different number of components can be discovered. As these cages tend to have common polyhedral shapes we can use this information to find linkers that can construct such shapes. 

### Development

The library is currently under development. As codebase gets bigger more documentation and examples will be provided.

[HD]: http://pubs.acs.org/doi/abs/10.1021/ic0202920
[Fig1]: https://github.com/kbsezginel/Moleidoscope/blob/master/doc/Figures/Fig1.PNG "Supramolecular Structure Generation"
[Fig2]: https://github.com/kbsezginel/Moleidoscope/blob/master/doc/Figures/Fig2.PNG "Supramolecular Cages"
