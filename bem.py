# %load_ext autoreload
# %autoreload 2
import sys
import logging, os
from time import time
import numpy as np
   
import matplotlib.pyplot as plt
from multiprocessing import Pool # rather than import multiprocessing to work with jupyter

sys.path.append('../../')   # add path of package "bem" to search list.
sys.path.append('../../../electrode/')   # add path of package "electrode" to search list.


from bem import Electrodes, Sphere, Mesh, Grid, Configuration, Result, Box
from bem.formats import stl

import numpy as np

scale = 40
prefix = 'htrapv8'

cmap = {
#     20083 : 'DC200', # base color?
#     32602 : 'GND',
# #     15366 : 'trench',
    31    : 'RF',
    25375 : 'DC21',
    32672 : 'DC1',
    32608 : 'DC2',
    32576 : 'DC3',
    32544 : 'DC4',
    32512 : 'DC5',
    32480 : 'DC6',
    32416 : 'DC7',
    32384 : 'DC8',
    32352 : 'DC9',
    32320 : 'DC10',
    5022  : 'DC11',
    5021  : 'DC12',
    5019  : 'DC13',
    5018  : 'DC14',
    5017  : 'DC15',
    5016  : 'DC16',
    5015  : 'DC17',
    5013  : 'DC18',
    5012  : 'DC19',
    5011  : 'DC20',
}

s = stl.read_stl(open(f'{prefix}.stl', "rb"))
stl.check_normals(*s[:2])
r = stl.stl_to_mesh(*s, scale=scale, rename=cmap) 
r = stl.split_by_normal(r)
mesh = Mesh.from_mesh(r)

# Plot mesh
fig, ax = plt.subplots(subplot_kw=dict(aspect="equal"), figsize=(16,10), dpi=300)
ax.set_xlabel("x/l",fontsize=10)
ax.set_ylabel("y/l",fontsize=10)
# ax.text(-1.5,7,"l = %d um"%(scale/1e-6),fontsize=12)
mesh.plot(ax)
# plt.savefig('fig.png', bbox_inches='tight')