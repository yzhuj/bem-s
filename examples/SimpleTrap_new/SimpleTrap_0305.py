#!/usr/bin/env python
# coding: utf-8
#here heref!
#more



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
import multiprocessing 
multiprocessing.set_start_method("fork")


import numpy as np

sys.path.append('../../')
from bem import Electrodes, Sphere, Mesh, Grid, Configuration, Result
from bem.formats import stl
import numpy as np
import copy
from helper_functions import *
from bemCol_lib.bemCol import rgb2stl, meshlab_bemCol, bemCol_dict, bemCol




# ### Import STL geometry file
# base file name for outputs and inputs is the script name

stl_file_in = "htrapf"
# scale to natural units (ion height)
# this seems not right to me- I feel like the ion-to-electrode distance is own for a spherical
# electrode geometry
leng = 1
factor = 1
# 72/leng
scale = 1e-3   # Distance from ion to electrode is 40 um.
use_stl = True

# we use a finite difference method to test converge:
# we set the mesh double finer, if the change of field is small, then we regard it converge
test_cvg = False


mesh,s_nta = load_file(Mesh,Electrodes,stl_file_in,scale,use_stl)
# The formal rename of electrode. Assign each electrode a string name instead of its color coding. Use the numbers you get above.
# `stl.stl_to_mesh()` prints normal vectors (different faces) in each electrode.

# print color with rename: 'bem1','bem2','uk1','uk2', etc
# s_nta[2] returns a set of attributes
print('test num color:',np.array(list(set(s_nta[2]))))
ele_col = bemCol(np.array(list(set(s_nta[2]))),('fusion360','export_stl'))
ele_col.print_colors_to_name()

ele_col.set_color_name(color = 'bem1',name = 'DC1')
ele_col.set_color_name(color = 'bem2',name = 'DC2')
ele_col.set_color_name(color = 'bem3',name = 'DC3')
ele_col.set_color_name(color = 'bem4',name = 'DC4')
ele_col.set_color_name(color = 'bem5',name = 'DC5')
ele_col.set_color_name(color = 'bem6',name = 'DC6')
ele_col.set_color_name(color = 'bem7',name = 'DC7')
ele_col.set_color_name(color = 'bem8',name = 'DC8')
ele_col.set_color_name(color = 'bem9',name = 'DC9')
ele_col.set_color_name(color = 'bem10',name = 'DC10')
ele_col.set_color_name(color = 'bem11',name = 'DC11')
ele_col.set_color_name(color = 'bem12',name = 'DC12')
ele_col.set_color_name(color = 'bem13',name = 'DC13')
ele_col.set_color_name(color = 'bem14',name = 'DC14')
ele_col.set_color_name(color = 'bem15',name = 'DC15')
ele_col.set_color_name(color = 'bem16',name = 'DC16')
ele_col.set_color_name(color = 'bem17',name = 'DC17')
ele_col.set_color_name(color = 'bem18',name = 'DC18')
ele_col.set_color_name(color = 'bem19',name = 'DC19')
ele_col.set_color_name(color = 'bem20',name = 'DC20')
ele_col.set_color_name(color = 'bem21',name = 'DC21')
ele_col.set_color_name(color = 'bem30',name = 'DC0')
ele_col.set_color_name(color = 'bem25',name = 'RF')

ele_col.print_drop_colors()

print(len(s_nta), type(s_nta),"\n")
# s_nta is a length 3 tuple. (normal, triangle, attribute)
# Normal direction of each triangle, three vetice of triangles, coding number of colors.
print("Triangles:",len(s_nta[0]),"\nColors:",len(s_nta[2]),"\n")    # This isn't right.

# stl_to_mesh() only assigns names and does scaling, doing no triangulation to stl mesh.
# "scale=scale/1e-6" only scales dimensionless scale/1e-6.    1e-6: if stl uses micron as unit.
color_dict = ele_col.result_dict
print('test_rename',color_dict)
mesh = Mesh.from_mesh(stl.stl_to_mesh(*s_nta, scale=1,
    rename=color_dict, quiet=False))



# ### Generate triangle mesh with constraints
#
# The meshes are 2-dimensional triangles on the surface of electrodes. The region enclosed by constraint shape can have finer mesh. Triangulation is done by `triangle` C library.
#there are all in units of mm now (Ben S. feb 2022)
xl = (3.1)*72*1e-3
yl = -0.051*72*1e-3
zl = 1.06*72*1e-3
rad = 5*72*1e-3
size = 100.0
inside=2e-4
outside=2e-3
# set .1 max area within 3
# areas_from_constraints specifies sphere with finer mesh inside it.
mpl.rcParams['lines.linewidth'] = 0.2 ###########

plot_mesh(xl,yl,mesh,scale,'fig1.png')

mesh.triangulate(opts="q10Q",new = False)
plot_mesh(xl,yl,mesh,scale,'fig2.png')
 # "inside", "outside" set different mesh densities.Q
# mesh.areas_from_constraints(Sphere(center=np.array([xl,yl,zl]),
#            radius=10*factor, inside=0.1*factor**2, outside=1000))
# mesh.triangulate(opts="",new = False)
print('--------------')
mesh.areas_from_constraints(Sphere(center=np.array([xl,yl,zl]),
           radius=rad, inside=inside, outside=outside))
# # retriangulate quality and quiet with areas
mesh.triangulate(opts="q20Q",new = False)
plot_mesh(xl,yl,mesh,scale,'fig3.png')
if test_cvg:
    mesh_fine = copy.deepcopy(mesh)
    mesh_fine.areas_from_constraints(Sphere(center=np.array([xl,yl,zl]),radius=rad, inside=inside/2, outside=outside/2))
    mesh_fine.triangulate(opts="q25Q",new = False)


# save base mesh to vtks
print("Output vtk:",os.path.abspath("./"+stl_file_in+".vtk"))    # output path



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
Lx, Ly, Lz = 0.100,0.100,0.100 # in the unit of scaled length l
sx, sy, sz = s, s, s

vtk_out = "vtks/htrapf"
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

jobs = list(Configuration.select(mesh,"DC.*","RF"))    # select() picks one electrode each time.

if test_cvg:
    jobs_fine = list(Configuration.select(mesh_fine,"DC.*","RF"))
else:
    jobs_fine = list(range(len(jobs)))

# run the different electrodes on the parallel pool
pmap = multiprocessing.Pool().map # parallel map
#pmap = map # serial map
t0 = time()

def run_map():
    list(pmap(run_job, ((jobs[i], grid, vtk_out, test_cvg , jobs_fine[i]) for i in range(len(jobs)))))
    print("Computing time: %f s"%(time()-t0))
    # run_job casts a word after finishing each electrode.

run_map()

fout = 'htrap_simulation_1_el3.5'
write_pickle(vtk_out,fout,grid)

# f = open('./' + 'gridExample' + '.pkl', 'wb')
# pickle.dump(grid, f, -1)
# f.close()


