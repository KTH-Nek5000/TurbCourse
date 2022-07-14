# Data-Driven Methods for Nek5000 Simulations
Python code for relevant data-driven techniques with I/O designed for Nek5000. So far, the available methods are:
  - Proper Orthogonal Decomposition (POD)
  - Dynamic Mode Decomposition (DMD), removed for this short tutorial

## General:
  - Nek5000 V19
  - Python 3 (anaconda)

## Usage:
  - Define the code settings in `input.txt`
  - Define the path to MODULES folder and input.txt file in `ddmMain.py`
  - The mass matrix has to be saved as `bm1casename0.f00001` and saved in the x-component
  - Please use this forked [pymech](https://github.com/danielemassaro/pymech) release

## Data:
  - The minimal-channel data can be found here on OneDrive [channel.tar](https://kth-my.sharepoint.com/:u:/g/personal/philipps_ug_kth_se/EbMQ_Pw2rQJMqez5KuK_eCABxYQ3zCC8LUZPB4kX-bvRgw?e=eP4lA1), size 1.5GB.
