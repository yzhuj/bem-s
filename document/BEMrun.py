import sys
sys.path.append('./helper_functions')
sys.path.append('./helper_functions/bemColors_lib')
sys.path.append('../../Electrodes')
from helper_functions import plot_mesh
from bemColors import bemColors
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from bem import Electrodes, Sphere, Mesh, Grid, Configuration, Result
from bem.formats import stl
from multipoles import MultipoleControl
from plottingfuncns import *
from time import time
from helper_functions import run_job, write_pickle
from BEMrun_config import *
from joblib import Parallel, delayed, cpu_count


# load stl file
from bem.formats import stl
s_nta = stl.read_stl(open(filename, "rb"))
print(f'LOADED {filename}')

# print colors needs to be named
bemcolors = bemColors(np.array(list(set(s_nta[2]))),('fusion360','export_stl'))
bemcolors.print_stl_colors()

for e in electrodes:
    bemcolors.color_electrode(*e)

# print colors still with no name. These meshes will be neglected in the following codes
bemcolors.drop_colors()

# read stl into mesh with electrode names
# unnamed meshes will not be imported at all
mesh = Mesh.from_mesh(stl.stl_to_mesh(*s_nta, scale=1, rename=bemcolors.electrode_colors, quiet=True))

mpl.rcParams['lines.linewidth'] = 0.2 

# plot what the stl looks like
if SHOW_PLOTS:
	plot_mesh(xl,yl,mesh,mesh_unit, title=filename, fout='fig1.png', save=SAVE_PLOTS, dpi=1000)

# Triangulation documentation: https://www.cs.cmu.edu/~quake/triangle.switch.html

# first triangulation 
mesh.triangulate(opts="q3Q",new = False)
if SHOW_PLOTS:
	plot_mesh(xl,yl,mesh,mesh_unit,title='first triangulation', fout='fig2.png', save=True, dpi=1000)

# second triangulation

# here we define a spherical constriant zone:
# units are in mesh_units
# inside and outside that sphere, triangle density is:
inside=10e-4
outside=10
rad = 10*75*1e-3


# areas_from_constraints specifies sphere with finer mesh inside it.
###############################################
# mesh.areas_from_constraints(Sphere(center=np.array([xl,yl,zl]),radius=rad, inside=inside, outside=outside))
# mesh.triangulate(opts="q30Q",new = False)
# if SHOW_PLOTS:
# 	plot_mesh(xl,yl,mesh,mesh_unit,title='final triangulation', fout='fig3.png', save=True, dpi=1000)

#proceed = input('Proceed with electrostatics? [Enter] to continue.')
proceed = 'n'
if proceed!='':
	print('EXITING')
	sys.exit()


# generate grid for BEM simulation
# grid to evalute potential and fields atCreate a grid in unit of scaled length mesh_unit. Only choose the interested region (trap center) to save time.
# I DON'T THINK n DOES ANYTHING HERE ######
n, s = 100, 0.015
#Lx, Ly, Lz = 2, 0.200 ,0.500 # in the unit of scaled length mesh_unit
Lx, Ly, Lz = .2, 0.200 ,0.200 # in the unit of scaled length mesh_unit
sx, sy, sz = s, s, s

prefix = "./pkls/" + filename.split('.')[0]
print('PREFIX: ', prefix)

# ni is grid point number, si is step size. Thus to fix size on i direction you need to fix ni*si.
nx, ny, nz = [2*np.ceil(L/2.0/s).astype('int') for L in (Lx, Ly, Lz)]
print("SIZE/l:", Lx, Ly, Lz)
print("STEP/l:", sx, sy, sz)
print("SHAPE (number of grid points):", nx, ny, nz)
grid = Grid(center=(xl,yl,zl), step=(sx, sy, sz), shape=(nx, ny, nz))
# Grid center (nx, ny ,nz)/2 is shifted to origin
print("GRID ORIGIN/l:", grid.get_origin())

# select() picks one electrode each time.
jobs = list(Configuration.select(mesh,"DC.*","RF"))    

t0 = time()

# without multiprocessing
#####################################################
if USE_MULTIPROCESSING:
    n_jobs = cpu_count()
    print(f'Running with {n_jobs} multiprocessing jobs')
    def run_map():
        Parallel(n_jobs=n_jobs)(delayed(run_job)((job, grid, prefix)) for job in jobs)
    run_map()
else:
    for job in jobs:
    	run_job((job, grid, prefix))

print(f"COMPUTING TIME: {time()-t0} s")

fout = open(prefix + '_grid.pkl', 'wb')
pickle.dump((mesh_unit,xl,yl,zl,grid), fout)
fout.close()
