Moleidoscope (Molecular Kaleidoscope)
=====================================

About
-----

<p align="center">
<img src="img/molecular-decomposition-animation.gif">
</p>

A key challenge in chemistry is to design molecules with given shapes and sizes comprising many components. One approach to designing structures with such a large number of components is to leverage the typically highly symmetrical nature of supramolecular structures (Figure 1). “Moleidoscope” interactively applies point group symmetry operations to generate hypothetical supramolecular structures in silico. Starting with a simple organic compound, the molecule is replicated using mirror operations, and by selecting multiple mirrors (3D planes), prisms can be formed to have multiple copies of the molecule oriented symmetrically. Moreover, as is the case with kaleidoscope, by rotating these mirrors even more copies can be obtained with many different symmetries.

---

<p align="center">
<img src="img/molecular-build-procedure.PNG">
</p>


**Figure 1:** Supramolecular structure generation process. First a library of supramolecular structures are decomposed into their building blocks which are simplified as building blocks. Using this library of building blocks new structures can be discovered by assembling these building blocks in different ways.

---

“Moleidoscope” performs these symmetry operations by selecting different organic compounds that can be found in [HostDesigner][HD], a software developed by Dr. Benjamin Hay in Oak Ridge National Laboratory. Also, it acts as a Python wrapper for HostDesigner which can be used to find linking fragments between the arrangement of molecules. In this way, different hypothetical supramolecular structures can be discovered.

### Supramolecular Cages

<p align="center"> <img src="img/supramolecular-cages.PNG"> </p>

Using _moleidoscope_ supramolecular cages with different number of components can be discovered. As these cages tend to have common polyhedral shapes we can use this information to find linkers that can construct such shapes.

### Setup

Moleidoscope requires Python ≥ 3.5.1 and [HostDesigner](https://github.com/kbsezginel/HostDesigner).

You can install Moleidoscope by cloning the repository and running setup.py as follows:

```bash
git clone https://github.com/kbsezginel/Moleidoscope.git
cd Moleidoscope
python setup.py install
```

The library is currently under development. As codebase gets bigger more documentation and examples will be provided.

---

[HD]: http://pubs.acs.org/doi/abs/10.1021/ic0202920
