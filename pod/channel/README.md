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

## Mirroring instructions
For example, let consider the flow around a 2D cylinder. The case is symmetric w.r.t. to x-axis. 
   1. In `$LOCALPATH/Nek5000/core/gfldr.f` modify the subroutine `gfldr(sourcefld)` as shown below. Before performing the interpolation, the read source mesh coordinates are swapped doing `ym1s=-ym1s`. In this way, we perform an interpolation exactly on the same mesh, but the GLL point at (x,y) is interpolated in (x,-y). Moreover, the order of GLL points and elements is maintained. This is fundamental when the snapshots matrix is assembled in the python code.
![`gfldr`](./docsrc/pics/gfldr_mirroring.png?style=centerme)
   2. Define in `case_name.usr` the usrchk as follows. We loop over all snapshots, performing the interpolation (as described above) and changing the sign of the v-velocity component. In this way, the mirroring operations are completed.
![`usrchk`](./docsrc/pics/usrchk_mirroring.png?style=centerme)
   3. Compile your code.
   4. Move the collected snapshots in the run-nek folder. Set in the `case_name.par` the parameter userParam01, which is the number of fields to mirror. Run your code and visualise the mirrored fields which have been generated.
