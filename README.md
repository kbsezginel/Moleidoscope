# Moleidoscope (Molecular Kaleidoscope)

## About

A key challenge in chemistry is to design molecules with given shapes and sizes comprising many components. One approach to designing structures with such a large number of components is to leverage the typically highly symmetrical nature of supramolecular structures (Figure 1). “Moleidoscope” interactively applies point group symmetry operations to generate hypothetical supramolecular structures in silico. Starting with a simple organic compound, the molecule is replicated using mirror operations, and by selecting multiple mirrors (3D planes), prisms can be formed to have multiple copies of the molecule oriented symmetrically. Moreover, as is the case with kaleidoscope, by rotating these mirrors even more copies can be obtained with many different symmetries. 
***
![alt text](https://github.com/kbsezginel/Moleidoscope/blob/master/doc/Figures/Fig1.PNG "Supramolecular Structure Generation")
**Figure 1:** Supramolecular structure generation process. First a library of supramolecular structures are decomposed into their building blocks which are simplified as building blocks. Using this library of building blocks new structures can be discovered by assembling these building blocks in different ways.
***
“Moleidoscope” performs these symmetry operations by selecting different organic compounds that can be found in HostDesigner,<sup>1</sup> a software developed by Dr. Benjamin Hay in Oak Ridge National Laboratory. Also, it acts as a Python wrapper for HostDesigner which can be used to find linking fragments between the arrangement of molecules. In this way, different hypothetical supramolecular structures can be discovered.

<sup>1</sup> Hay, B. P., & Firman, T. K. (2002). HostDesigner: A program for the de novo structure-based design of molecular receptors with binding sites that complement metal ion guests. Inorganic chemistry, 41(21), 5502-5512.
