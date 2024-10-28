# XTIgA-Waves
**Add-On to GeoPDEs for Stable Space-Time Isogeometric Wave Equation**

This repository contains the experimental source code to reproduce the stabilized space-time isogeometric method for the acoustic wave equation of
[An unconditionally stable space–time isogeometric method for the acoustic wave equation](https://arxiv.org/abs/2303.07268),
and the unconditionally stable space-time isogometric method of 
[Unconditionally stable space-time isogeometric discretization for the wave equation in Hamiltonian formulation](https://arxiv.org/abs/2303.07268).

## Installing the package
```bash
git clone --recursive https://github.com/XTIgA-Waves/XTIgA-Waves
```

### Dependencies

**MATLAB installation of GeoPDEs**

  * Download and uncompress the file [GeoPDEs_full.tar.gz](https://rafavzqz.github.io/geopdes/download/). 
  * Uncompress and untar the NURBS and the GeoPDEs packages, in the files
    nurbs-<version>.tar.gz, geopdes-<version>.tar.gz, respectively.
  * Add the generated folders to the path, including their subfolders. 
    You can do this by typing in in the command window:

    addpath (genpath ('nurbs'));
    
    addpath (genpath ('geopdes'));
    
    addpath (genpath ('geopdes_hierarchical'));
    
  * Install the mex-files for the "nurbs" package (OPTIONAL):
     * uncompress and untar the file nurbs_mex_files.tar.gz in the folder "nurbs/inst" of the nurbs package;
     * go to the folder "nurbs/inst" and run the script file 'compile'.
       
     This will compile the files and save the nurbs package to your MATLAB path.

## News
Nov, 2024: added code related to the numerical tests and verifications of [Unconditionally stable space-time isogeometric discretization for the wave equation in Hamiltonian formulation](https://arxiv.org/abs/2303.07268).
*Nov, 2023: First release of `XTIgA-Waves`. 

## Papers using the code
* *An unconditionally stable space–time isogeometric method for the acoustic wave equation*  
Sara Fraschini, Gabriele Loli, Andrea Moiola, Giancarlo Sangalli  
[![arXiv](https://img.shields.io/badge/arXiv-2303.07268-b31b1b.svg)](https://arxiv.org/abs/2303.07268)

* *Unconditionally stable space-time isogeometric discretization for the wave equation in Hamiltonian formulation*  
Matteo Ferrari, Sara Fraschini, Gabriele Loli, Ilaria Perugia
[![arXiv](https://img.shields.io/badge/arXiv-2303.07268-b31b1b.svg)](https://arxiv.org/abs/2303.07268)


If you are using `XTIgA-Waves` in your academic work, please consider citing 
```
Fraschini, S., Loli, G., Moiola, A., & Sangalli, G. (2023).
An unconditionally stable space-time isogeometric method for the acoustic wave equation.
arXiv preprint arXiv:2303.07268
```
```
Ferrari, M. Fraschini, S., Loli, G., Perugia, I. (2024).
Unconditionally stable space-time isogeometric discretization for the wave equation in Hamiltonian formulation
arXiv preprint arXiv:2303.07268
```
