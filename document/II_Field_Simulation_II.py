from joblib import Parallel, delayed, cpu_count

from time import time
import pickle
import numpy as np
from bem import Electrodes, Sphere, Mesh, Grid, Configuration, Result
import multiprocessing 
import sys


sys.path.append('./helper_functions')
from helper_functions import run_job, write_pickle

fout_name = './intermediate_results/testmesh.pkl'
with open(fout_name,'rb') as f:
    mesh_unit,xl,yl,zl,mesh = pickle.load(f) # import results from mesh processing

# grid to evalute potential and fields atCreate a grid in unit of scaled length mesh_unit. Only choose the interested region (trap center) to save time.
n, s = 100, 0.100
Lx, Ly, Lz = 0.100,0.100,0.100 # in the unit of scaled length mesh_unit
sx, sy, sz = s, s, s

vtk_out = "./intermediate_results/.vtks/htrapf"
print("done")

# ni is grid point number, si is step size. Thus to fix size on i direction you need to fix ni*si.
nx, ny, nz = [2*np.ceil(L/2.0/s).astype('int') for L in (Lx, Ly, Lz)]
print("Size/l:", Lx, Ly, Lz)
print("Step/l:", sx, sy, sz)
print("Shape (grid point numbers):", nx, ny, nz)
grid = Grid(center=(xl,yl,zl), step=(sx, sy, sz), shape=(nx, ny, nz))
# Grid center (nx, ny ,nz)/2 is shifted to origin
print("Grid origin/l:", grid.get_origin())


# ## (2) run jobs
# 
# evaluate electric potential for each configurations: one electrode at 1V, the rest 0V. Set `pmap` as `multiprocessing.Pool().map` for parallel computing. For serial map, use `map`.

# In[5]:


jobs = list(Configuration.select(mesh,"DC.*","RF"))    # select() picks one electrode each time.


# run the different electrodes on the parallel pool
# pmap = multiprocessing.Pool().map # parallel map
#pmap = map # serial map
t0 = time()

# using multiprocessing
#def run_map():
#    list(pmap(run_job, ((job, grid, vtk_out) for job in jobs)))
#    print("Computing time: %f s"%(time()-t0))
    # run_job casts a word after finishing each electrode.

# without multiprocessing
for job in jobs:
	run_job((job, grid, vtk_out))
	# print(job.name)

# multiprocessing with joblib
# n_jobs = cpu_count()
# print(f'Running with {n_jobs} multiprocessing jobs')
# def run_map():
    # Parallel(n_jobs=n_jobs)(delayed(run_job)((job, grid, vtk_out)) for job in jobs)
# run_map()

print("Computing time: %f s"%(time()-t0))
    
#multiprocessing.set_start_method("fork", force=True)

fout = './intermediate_results/test_field_result'
write_pickle(vtk_out,fout,grid)

print('Done')





