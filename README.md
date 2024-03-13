# MD_and_NMR_Relaxometry

## Description

This repo contains the simulation setups and analysis script for the publication on ubiquitin by Champion et al.
All simulations were performed with GROMACS and the topology files and simulation setups can be found in the folder `simulations`. No trajectories are provided in this repo but can be provided upon request. The simulations are grouped in folders with the names of the different force field parameters used. For further information on the simulation setups, please refer to the publication.

**Unraveling motion in proteins by combining NMR relaxometry and molecular dynamics simulations: A case study on ubiquitin**  
Candide Champion, Marc Lehner, Albert A. Smith, Fabien Ferrage, Nicolas Bolik-Coulon, Sereina Riniker   
J. Chem. Phys. 160, 104105 (2024)  
https://doi.org/10.1063/5.0188416

## Installation

The analysis scripts are all written in Python. A conda environment file is provided as `environment.yml`. To install the environment, run the following command in the terminal:

```bash
conda env create -f environment.yml
```

This will create a conda environment called `MD_and_NMR_Relaxometry`. To activate the environment, run:

```bash
conda activate MD_and_NMR_Relaxometry
```

Additionally, the python package pyDR (by A. Smith) is required, which can be found [here](https://github.com/alsinmr/pyDR). The package can be installed by running:

```bash
pip install git+https://github.com/alsinmr/pyDR.git
```

The python package Energy based Clustering can be found here: [here](https://github.com/rinikerlab/EnergyBasedClustering). The package can be installed with:

```bash
pip install git+https://github.com/rinikerlab/EnergyBasedClustering.git
```

## Usage

The analysis scripts are all written in Python, mostly as Jupyter notebooks. The notebooks are located in the folder `analysis`. The notbooks can be used to reproduce the figures in the publication.
