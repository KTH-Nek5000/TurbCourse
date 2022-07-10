## A guide to apply POD to the snapshot of a 2d mixing layer simulated by Nek5000.


### Dependencies: 
The POD code is written in `Python` and requires the following libraries to be installed:
    * [`numpy`](https://numpy.org/)
    * [`matplotlib`](https://matplotlib.org/)
    * [`pickle`](https://docs.python.org/3/library/pickle.html#module-pickle)    
    * [`pymech`](https://pymech.readthedocs.io/en/latest/) (optional): to read the Nek5000 fields and create a `pickle` database. 

To run the notebook, a set of snapshot data of a 2D mixing layer simulated by `Nek5000` are provided in folder `./mixlay2d_snaps/`. 
These snapshots are `mixlay0.f*`. To read them, [`pymech`](https://pymech.readthedocs.io/en/latest/) is needed to be installed. 
As an alternative, a [`pickle`](https://docs.python.org/3/library/pickle.html#module-pickle) databse is created from the same set of snapshots which can be readily be used without having `pymech` installed. 

