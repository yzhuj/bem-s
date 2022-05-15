# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## 2022-04 Wenhao He [Unreleased]

### Fixed
- errors in triangulate wrap functions
- For 2D trap, parallel planes no longer classified as same plane
- combine same points in Mesh object, thus eliminate sharp triangles
- make triangles counter clock wise for triangulation package
- errors on CentOS(HPC)

### Change
- abandon mayavi for much more convinient way of maitaining the code
- new way to name the electrode
- add a set of standard colors
- add a document
- combine BEM and Electrode togethor


### to do
- Simulating different trap geometries




## 2022-04 Wenhao He [Unreleased]

### Fixed
- change files in fastlap folder same as master branch to avoid declaration error
- In SimpleTrap_0305.py, add a setting on multiprocess because of its changed realization method in Python 3.9.12
- veryfied the unit is mm when `solve singulatrities` and `simulate` manually.

### Change
- simplify the environment set up process(user friendly), with one command"pip install ."
- tested compatible with latest vtk==9.0.3, latest mayavi==4.7.4, and old PyQt5==5.15.6 on python==3.9.12 on both Linux(Ubuntu 20.04.3) and macOS(Monterey 12.1 intel)
- In SimpleTrap_0305.py, add testing convergence(`test_cvg`). Making mesh area bound half and compute relative error of potential. 
- seperate readme into readme + changelog + COPYING


### to do
- I think we'd better break away from package mayavi(tvtk.api) to avoid its problem of compatibility with new version python. There are two benefits: 1. save a lot of effort for maintaining our code 2.If running our code in super computer(probably CentOs), being without mayavi would be more friendly. Now our code invoke package mayavi(tvtk.api) for two reasons:1. for data visulization 2. processing vtk file. For data visulization, we can use other packages, the real challenge for abandoning mayavi is to solve the second problem. There are two probable solution:1.use an alternative for processing vtk files 2. abandon .vtk files also, and generate .pkl file directly.








## 2019-03 Wance Wang [Unreleased]

### Change
- setup process for python 3.6:

Create a python environment "ele36" for bem package:
```
    $ conda create -n ele36 python=3.6  
    $ conda install -n ele36 jupyter scipy matplotlib cython cvxopt apptools envisage  
    $ source activate ele36
    $ pip install mayavi
```
Or you can reproduce my environment by using file `ele36pip_2019.yml`.
```
    $ conda env create -f ele36pip_2019.yml  
    $ source activate ele36  
```
For building `fastlap` and `triangle` C extensions:
```
    $ python setup.py build_ext --inplace
```
(When you want to rebuild it, add a `--force` option.)




