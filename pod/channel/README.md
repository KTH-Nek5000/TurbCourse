# Data-Driven Methods for Nek5000 Simulations
Python code for relevant data-driven techniques with I/O designed for Nek5000. So far, the available methods are:
  - Proper Orthogonal Decomposition (POD)
  - Dynamic Mode Decomposition (DMD)

## General:
  - Nek5000 V19
  - Python 3

## Usage:
  - Define the code settings in `input.txt`
  - Define the path to MODULES folder and input.txt file in `ddmMain.py`
  - The mass matrix has to be saved as `bm1casename0.f00001` and saved in the x-component
  - Please use this forked [pymech](https://github.com/danielemassaro/pymech) release
