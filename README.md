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

Nov, 2023: First release of `XTIgA-Waves`. 

## Papers using the code
* *An unconditionally stable space–time isogeometric method for the acoustic wave equation*  
Sara Fraschini, Gabriele Loli, Andrea Moiola, Giancarlo Sangalli  
[<img class="publication-brand-image" src="https://sdfestaticassets-eu-west-1.sciencedirectassets.com/prod/4b190b964b2415d63e7a2050ee5b17f5f8cbb4eb/image/elsevier-non-solus.png" alt="Elsevier" width="50">](https://www.sciencedirect.com/science/article/pii/S0898122124002773)


* *Space–time isogeometric analysis: a review with application to wave propagation*  
Gabriele Loli, Giancarlo Sangalli
[<img class="publication-brand-image" src="https://media.springernature.com/w316/springer-static/cover-hires/journal/40324?as=webp" alt="SeMA" width="50">](https://doi.org/10.1007/s40324-025-00391-x)




* *Unconditionally stable space-time isogeometric discretization for the wave equation in Hamiltonian formulation*  
Matteo Ferrari, Sara Fraschini, Gabriele Loli, Ilaria Perugia  
[![arXiv](https://img.shields.io/badge/arXiv-2411.00650-b31b1b.svg)](https://arxiv.org/abs/2411.00650)


If you are using `XTIgA-Waves` in your academic work, please consider citing 
```
Fraschini, S., Loli, G., Moiola, A., & Sangalli, G. (2024).
An unconditionally stable space-time isogeometric method for the acoustic wave equation.
Computers & Mathematics with Applications, Volume 169, 2024, Pages 205-222, ISSN 0898-1221.
https://doi.org/10.1016/j.camwa.2024.06.009.
```
```
Loli, G., Sangalli, G. (2025).
Space–time isogeometric analysis: a review with application to wave propagation
SeMA, 2025.
https://doi.org/10.1007/s40324-025-00391-x
```
```
Ferrari, M., Fraschini, S., Loli, G., Perugia, I. (2024).
Unconditionally stable space-time isogeometric discretization for the wave equation in Hamiltonian formulation
arXiv preprint arXiv:2411.00650
```

# Contributions 
Lead Developers - *Sara Fraschini* and *Gabriele Loli*, Core development and code management 

Developers  - *Matteo Ferrari*, Development of the folder `verifications`

Idea Contributors - *Andrea Moiola, Ilaria Perugia* and *Giancarlo Sangalli*
