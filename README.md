# Description
IDECAMB is a CAMB and CosmoMC patch for interacting dark energy (IDE) models. It provides a unified interface for the widely studied coupled quintessence (CQ) and the coupled fluid (CF) models. The perturbation evolutions of the CF models are handled by the extended parametrized post-Friedmann (PPF) approach to avoid the possible large-scale instability.

# Installation 
A simple way to install this patch is just copying all the files into the [CosmoMC-planck2018 package](https://github.com/cmbant/CosmoMC/tree/planck2018) and overwriting the same files. Or you can modify the CosmoMC package with reference to the file modification_manual.f90

# Usage and running
Modify test_ide.ini, then run ./cosmomc test_ide.ini

# How to cite us
If you use IDECAMB, please cite its pre-print, [arXiv:1404.5220](https://arxiv.org/abs/1404.5220) and [arXiv:2306.01593](https://arxiv.org/abs/2306.01593).
