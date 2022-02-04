#!/usr/bin/env python
# coding: utf-8
#here heref!
#more
# In[1]:


#   bem: triangulation and fmm/bem electrostatics tools
#
#   Copyright (C) 2011-2012 Robert Jordens <jordens@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
from time import time
import numpy as np
import resource
import matplotlib as mpl
import matplotlib.pyplot as plt
from multiprocessing import Pool

sys.path.append('/Users/Ben/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/ionLifetimes/bem')
sys.path.append('/Users/Ben/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/ionLifetimes/bem/examples')
sys.path.append('/Users/Ben/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/ionLifetimes/electrode')
import numpy as np
from helper_functions import *
from bem import Electrodes, Sphere, Mesh, Grid, Configuration, Result
from bem.formats import stl
from trap_library import *
import numpy as np


# ### Import STL geometry file
# base file name for outputs and inputs is the script name

prefix = "htrap13-14_6-11-12-4-5-8-gnd"
suffix = ""
# scale to natural units (ion height)
# this seems not right to me- I feel like the ion-to-electrode distance is own for a spherical
# electrode geometry
leng = 1
factor = 1
# 72/leng
scale = 1e-3   # Distance from ion to electrode is 40 um.
use_stl = True

mesh,s_nta = load_file(Mesh,Electrodes,prefix,scale,use_stl)
# The formal rename of electrode. Assign each electrode a string name instead of its color coding. Use the numbers you get above.
# `stl.stl_to_mesh()` prints normal vectors (different faces) in each electrode.


print(len(s_nta), type(s_nta),"\n")
# s_nta is a length 3 tuple. (normal, triangle, attribute)
# Normal direction of each triangle, three vetice of triangles, coding number of colors.
print("Triangles:",len(s_nta[0]),"\nColors:",len(s_nta[2]),"\n")    # This isn't right.

# stl_to_mesh() only assigns names and does scaling, doing no triangulation to stl mesh.
# "scale=scale/1e-6" only scales dimensionless scale/1e-6.    1e-6: if stl uses micron as unit.

mesh = Mesh.from_mesh(stl.stl_to_mesh(*s_nta, scale=1,
    rename=el_colordict[prefix], quiet=False))

# ### Generate triangle mesh with constraints
#
# The meshes are 2-dimensional triangles on the surface of electrodes. The region enclosed by constraint shape can have finer mesh. Triangulation is done by `triangle` C library.

xl = 3.7*72*1e-3
yl = -0.051*72*1e-3
zl = 1.06*72*1e-3
# set .1 max area within 3
# areas_from_constraints specifies sphere with finer mesh inside it.
mpl.rcParams['lines.linewidth'] = 0.2
rad = 5*72*1e-3
size = 100.0
file_name = "el2(4-5-6-8-11-12-gnd_13-14).txt"
print(file_name)
 # "inside", "outside" set different mesh densities.
# mesh.areas_from_constraints(Sphere(center=np.array([xl,yl,zl]),
#            radius=10*factor, inside=0.1*factor**2, outside=1000))
# mesh.triangulate(opts="",new = False)
mesh.areas_from_constraints(Sphere(center=np.array([xl,yl,zl]),
           radius=rad, inside=2e-4, outside=1e-3))
# # retriangulate quality and quiet with areas
mesh.triangulate(opts="",new = False)
# save base mesh to vtks
# mesh.to_vtk(prefix+suffix)
# mesh.to_vtk(prefix+suffix)
print("Output vtk:",os.path.abspath("./"+prefix+suffix+".vtk"))    # output path

plot_mesh(xl,yl,mesh,scale)


# ### Main boundary element calculations
#
# In `run_job` function, `job` is `Configuration` instance and `grid` is discretirized spatial grid (not the mesh). The general workflow (also the routine of BEM method) are:
# 1. `solve_singularities()` solves charge distributions by iterative methods to make it consistent with one electrode at 1V and others at 0V (unit potentials). `adapt_mesh()` refines meshes adaptively to achieve certain precision while solving sigulartities.
# 2. Compute potentials on given grid points by `simulate()`, based on the charge distributions gotten previously.
# 3. Potential data of each unit potential are saved seperately to a `Result` instance, and also export to VTK files.
# 4. Return total accumulated charge per electrode in the end.
#
# Major calculations calls `fastlap` C library which uses a pre-conditioned, adaptive, multipole-accelerated algorithm for solving Laplace problem. Two parameters control multipole acceleration.
# + num_mom, the number of multipole
# + num_lev, the number of levels in the hierarchical spatial decomposition.
# num_lev=1 means direct computation without multipole acceleration. See fastlap ug.pdf and README.rst.

# Create a grid in unit of scaled length `l`. Only choose the interested region (trap center) to save time.
# For reference, to compute Seidelin trap, grid shape = (60, 60, 60) takes 266 s, while shape = (150, 150, 150) takes 3369 s.

# grid to evalute potential and fields atCreate a grid in unit of scaled length l. Only choose the interested region (trap center) to save time.
n, s = 100, 0.002
Lx, Ly, Lz = 0.150,0.075,0.075 # in the unit of scaled length l
sx, sy, sz = s, s, s

prefix = "htrapF_mega_short"+str(s)+"_size"+str(size)+""
# prefix = "htrap_13-14_6-gnd_11-12-gnd"

# os.mkdir(prefix)
suffix = ""
print("done")

# ni is grid point number, si is step size. Thus to fix size on i direction you need to fix ni*si.
nx, ny, nz = [2*np.ceil(L/2.0/s).astype('int') for L in (Lx, Ly, Lz)]
print("Size/l:", Lx, Ly, Lz)
print("Step/l:", sx, sy, sz)
print("Shape (grid point numbers):", nx, ny, nz)
grid = Grid(center=(xl,yl,zl), step=(sx, sy, sz), shape=(nx, ny, nz))
# Grid center (nx, ny ,nz)/2 is shifted to origin
print("Grid origin/l:", grid.get_origin())

# Calculation. Parallel computation `Pool().map`
# generate electrode potential configurations to simulate
# use regexps to match electrode names
jobs = list(Configuration.select(mesh, "DC.*","RF"))    # select() picks one electrode each time.
# run the different electrodes on the parallel pool
pmap = Pool().map # parallel map
# pmap = map # serial map
t0 = time()
def run_map():
    list(pmap(run_job, ((job, grid, prefix+suffix) for job in jobs)))
    print("Computing time: %f s"%(time()-t0))
    # run_job casts a word after finishing each electrode.
#
# run_map()
#s
#