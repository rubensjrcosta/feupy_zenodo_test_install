[![gammapy](https://img.shields.io/badge/powered%20by-gammapy-orange.svg?style=flat)](https://gammapy.org/)

# FeuPy

**FeuPy** is a Python application for very high-energy (VHE) gamma-ray analysis using the [Gammapy](https://docs.gammapy.org/1.3/) library. It facilitates the search for possible γ-ray counterparts to target sources and performs spectral model fitting across multiple wavelengths. **FeuPy** is designed with a focus on Cherenkov Telescope Array (CTA) users, providing tools for sensitivity analysis, observation simulations, and the computation of non-thermal radiation from relativistic particle populations.

## Features

- **Multi-Wavelength Spectral Model Fitting**: Perform model fitting across different energy bands using the Gammapy library.
  
- **Catalog Integration**: Incorporates data from multiple external catalogs:
  - [PSRCAT Pulsar Catalog](https://www.atnf.csiro.au/research/pulsar/psrcat/)
  - [VERITAS Catalog](https://github.com/VERITAS-Observatory/VERITAS-VTSCat/tree/main)
  - Additional observations from dedicated publications

- **Cherenkov Telescope Array (CTA) Sensitivity Analysis**: Provides built-in tools for conducting sensitivity analysis based on CTA performance. 
  - [CTA Sensitivity Analysis Tutorial](https://docs.gammapy.org/1.1/tutorials/analysis-1d/cta_sensitivity.html)
  
- **Observation Simulations**: Simulate observations and explore expected results under different scenarios. 
  - [Spectrum Simulation Tutorial](https://docs.gammapy.org/1.1/tutorials/analysis-1d/spectrum_simulation.html)
  
- **Non-Thermal Radiation Computation**: Supports the computation of non-thermal emission from relativistic particle populations using the [Naima](https://naima.readthedocs.io/en/latest/examples.html) library.

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/feupy.git
   cd feupy
2. The last conda commands will define the environment variable within the conda environment. Conversely, you might want to define the $PYTHONPATH environment variable directly in your shell with:
   export PYTHONPATH=/your-feupy-path/feupy
