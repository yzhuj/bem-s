# BEM
Wraps the fastlap BEM (boundary element method aka BIE boundary integral
equations) + FMM (fast multipole method) + Laplace kernel implementation
from [http://www.rle.mit.edu/cpg/research_codes.htm](http://www.rle.mit.edu/cpg/research_codes.htm). The original code is
in /doc/fastlap_fl-2.0-22oct96 for reference.

## Installation

For users:

Python version <= 3.9 is required. 

Anaconda environment management system is recommended:

    $ conda create -n ele39 python=3.9
    
    $ conda activate ele39

To set up (in folder of setup.py, setup.cfg, pyproject.toml and MANIFEST.in):

    $ pip install .

to use mayavi dependent packages(optional):

    $ pip install 'vtk<9.1' 'PyQt5<5.16'; pip install 'mayavi<4.8'

The workflow of our package is explained in examples/SimpleTrap/SimpleTrap_0305.py

5/16/2022: tested compatible with latest vtk==9.0.3, latest mayavi==4.7.4, and old PyQt5==5.15.6 on python==3.9.12 on Linux(Ubuntu 20.04.3) and macOS(Monterey 12.1 intel) and CentOS(without mayavi package). More tests are needed for Windows.

## Documents and tutorial
Document for the latest version can be found at ./examples/document, which explains our workflow in detail.

Codes in ./examples/SimpleTrap is also helpful.


## Notes

1. The package can read STL files. Please color different electrodes with different colors. This can be done in Inventor. When exporting to STL in Inventor, be sure to choose microns for the Units. The code scales the trap as if it is in units if microns. There are no curved surfaces, so Surface and Normal deviations are not important. Set max edge length to something reasonable or else triangles are too small or large (40 microns, the ion height from the surface, is a good choice). Choose a low aspect ratio and make sure to export the colors.
2. The meshing is now done preferably with Jonathan Shewchuk's triangle code. It is integrated into the Python package and compiled into a python module, see /triangle and /bem/pytriangle.pyx. After read STL file and convert it to mesh, further meshing can be done using `Mesh.triangulate()`. Please check 'opts' being passed to `Mesh.triangulate()`. It follows the same command line switch syntax as the original code. Frequently used are angle and size constraints: [https://www.cs.cmu.edu/~quake/triangle.quality.html](https://www.cs.cmu.edu/~quake/triangle.quality.html).
3. The meshing can be adaptive: for a given potential configuration, the charge distribution for a small initial mesh is solved and the mesh is refined such that each triangle in the final mesh contributes equally to the field at the observation point (since the triangles are also constrained by the geometry, min angle etc, the final number of triangles is typically larger than the desired one).
3. The meshing can also be done in Inventor, prior to using bem package. Please make sure to run `Mesh.gather()` if not using `Mesh.triangulate()`. 
3. The Fastlap code is also integrated and compiled into a python module (see /fastlap and /bem/fastlap.pyx and called directly without exporting and importing intermediate file formats.
4. Triangle areas can be constrained via differently shaped constraint fields (Box, Sphere, Cylinder or BYO).
5. Please read jupyter notebooks in /examples folder for more instructions.


## Changlog
Changes of historical versions are documented in CHANGELOG.md
