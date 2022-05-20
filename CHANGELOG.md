# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## 2022-05 Wenhao He

### Change
- simulate overhang htrap(3D electrode): original code(especially triangulation part) behaves bad for 3D trap, so a lot of changes are made in Bem.
    - unify same points in mesh structure
    - counter clockwise triangle when passed into the package
    - originally, parallel planes will be classified as one plane, which leads to error for 3D cases. Fix now
    - now we can combine triangles first and cut it, so sharp triangles are eliminated
    - fix the error: exceed float precision .......
- a more convenient way to load color (explained in the document)
    - a set of standard color
    - more regular color names
- simplify the environment set up process, with one command: pip install .
    - setup.py, setup.cfg, pyproject.toml and MANIFEST.in
- compatible with python 3.9
    - solve the problem of multiprocessing
    - solve the environment for mayavi PyQt vtk
    - get around mayavi in our workflow
        - our code can run on HPC(CentOs)
        - saves a lot of effort to maintain the code in the future
- others
    - make unit problem clear(micronmeter or millimeter are both ok now)
    - split readme into readme+changelog+copying
    - bugs in ./fastlap
    - other small bugs

- add a document: ./examples/document











