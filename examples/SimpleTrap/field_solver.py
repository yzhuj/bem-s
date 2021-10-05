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

from helper_functions import *
from bem import Electrodes, Sphere, Mesh, Grid, Configuration, Result
from bem.formats import stl
from trap_library import *
from SimpleTrap_0305 import *
# Contour plot of potential/pseudo-potential in 3 directions
# isocontour plot of RF pseudopotential radially from x (axial) direction

plot_RF(Result, prefix, suffix, grid)

# isocontour plot of DC potential from x (axial) direction
strs = "DC1 DC2 DC3 DC4 DC5 DC6 DC7 DC8 DC9 DC10 DC11 DC12 DC13 DC14 DC15 DC16 DC17 DC18 DC19 DC20 DC21".split()
plot_DC(Result, prefix, suffix, grid, strs, dir='x')

# isocontour plot of electrode potential (electrode profile) from z direction


# ### 3D plot of mesh and potential isocontour
# By mayavi GUI (seems to have problem now.)
# explore it in fancy 3D
# fire up a mayavi2 window showing base mesh, charges on final mesh
# and isosurfaces of the pseudopotential
Result.view(prefix + suffix, "RF")
# need to start the full eventloop for the window.
# close it to return control to the notebook
from pyface.api import GUI

GUI().start_event_loop()

# electrode is an another package in nist-ionstorage github. wwc
from electrode import System, GridElectrode, utils
import csv
import scipy.constants as ct

s = System()
# load the electrostatics results into a electrode.System()
# "DC0 DC1 DC2 DC3 DC4 DC5 DC6 DC7 DC8 DC9 DC10 DC11 DC12 DC13 DC14 DC15 DC16 DC17 DC18 DC19 DC20 DC21 RF
# strs = "DC1 DC2 DC3 DC4 DC5 DC6 DC7 DC8 DC9 DC10 DC11 DC12 DC13 DC14 DC15 DC16 DC17 DC18 DC19 DC20 DC21".split()
strs = "DC1 DC2 DC3 DC4 DC5 DC6 DC7 DC8 DC9 DC10 DC11 DC12 DC13 DC14 DC15 DC16 DC17 DC18 DC19 DC20 DC21".split()
excl = {
    "DC6": ["Null", 0],
    "DC14": [13, 12]
}
exp = 4
offset = 1

npl = 9 - offset
lamb = np.zeros((len(strs), npl))
i = 0

r = Result.from_vtk(prefix, "DC1")
e = GridElectrode.from_result(r, maxderiv=exp)
print("here")
print(np.shape(e.data[1]))
for name in strs:
    if name in excl:
        print(name)
        arlo = np.zeros(npl)
    else:
        r = Result.from_vtk(prefix, name)
        e = GridElectrode.from_result(r, maxderiv=exp)
        e.name = name
        s.append(e)
        arlo = np.array([])
        for p in np.arange(offset, exp):
            if name == 'RF':
                print("trip")
            else:
                sx = 2
                sy = 2
                sz = 2
                #             print(vx)
                pot = np.zeros(2 * p + 1)
                pt = np.array([xl, yl, zl])
            outvals = utils.cartesian_to_spherical_harmonics(np.transpose(e.data[p][nx // 2 - 1:nx // 2 + 1,
                                                                          ny // 2 - 1:ny // 2 + 1,
                                                                          nz // 2 - 1:nz // 2 + 1, :]))
            outvals = np.sum(np.sum(np.sum(outvals, 1), 1), 1)
            arlo = np.append(arlo, outvals)
    lamb[i] = arlo[0:npl]
    i = i + 1
print(np.shape(lamb))
l = 72e-6  # length scale
u = 100.  # peak rf voltage
o = 36.e6 * 2 * np.pi  # rf frequency
m = 40 * ct.atomic_mass  # ion mass
q = 1 * ct.elementary_charge  # ion charge
rf_scale = s.rf_scale(m, q, l, o)
r = Result.from_vtk(prefix + suffix, "RF")
p = r.pseudo_potential
# coefficient of pseudo-potential. See blakestad2010 Eq.(5.2). Run rf_scale() before other calculations.
e = GridElectrode.from_result(r)
e.name = "RF"
s.append(e)
s["RF"].rf = u  # peak rf voltage

#     print(outvals[0][2])
#     print(utils.cartesian_to_spherical_harmonics(e.data[0:1][5][5][5]))
#     e = GridElectrode.from_vtk("%s%s_"%(prefix,suffix)+name+".vtk",maxderiv=4)
# print(np.dot(np.transpose(lamb),np.transpose(np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]))))
# print(arlo)


# In[15]:


# u2old = np.array([-0.0178185, 0.0153893, 0.0434955, 0.473985, 3.0024, -3.33158, \
# -4.52325, -0.919275, -0.716663, -0.784725, -0.0289778, 0.0169508, \
# 0.0820125, 0.601035, 2.83178, -3.5178, -6.54885, -0.756, 0.840975, \
# -0.0394898, -0.90165])/7
u2old = np.array(
    [-0.054593, 0.057666, -0.66996, 0.05951, -0.057175, -0.058409, -0.047644, -0.038418, -0.031, -0.02519, -0.045262,
     -0.077023, -0.68039, -0.074042, -0.034745, -0.039145, -0.03209, -0.025515, -0.019994, -0.01579, -0.086705])
lambT = np.transpose(lamb)
# commando = np.linalg.inv(lambT)
commandoT = np.zeros((8, len(strs)))
u2vec = np.zeros(8)
u2vec[6 - offset] = 1
for i in np.arange(0, 8):
    B = np.zeros(npl)
    B[i] = 1
    A = np.linalg.lstsq(lambT, B, rcond=None)
    commandoT[i] = A[0] * (np.dot(lambT, u2old)[5])
lambT = lambT[0:8]
commando = np.transpose(commandoT)
from tabulate import tabulate

print("checking construction")
print(np.dot(lambT, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]))
print(arlo)

print("checking inverse")
print(tabulate(np.around(np.dot(commando, lambT))))

print("u2 from pinv")
u2fromInv = np.around(np.dot(commando, u2vec), 1)
# print(u2fromInv)
for stri in excl:
    print(stri)
    idx = excl[stri][0]
    new = excl[stri][1]
    if idx != "Null":
        commandoT[:, idx] = commandoT[:, new]
l1 = commandoT.flatten()
lmid = l1
lmid[63:84] = l1[147:168]
lmid[84:105] = l1[105:126]
lmid[105:126] = l1[63:84]
lmid[126:147] = l1[84:105]
lmid[147:168] = l1[126:147]
# lmid = lamb.flatten()
print("u2 from indexing")
u2 = np.around(commando[:, 5], 5)
# el3 solution
u2 = np.array(
    [-0.054593, 0.057666, -0.66996, 0.05951, -0.057175, -0.058409, -0.047644, -0.038418, -0.031, -0.02519, -0.045262,
     -0.077023, -0.68039, -0.074042, -0.034745, -0.039145, -0.03209, -0.025515, -0.019994, -0.01579, -0.086705])
# el7 solution
print("getting u2 again")
print(np.around(np.dot(lambT, u2) / (np.dot(lambT, u2old)[5]), 3))
print(np.linalg.norm((np.dot(lambT, u2) - [0, 0, 0, 0, 0, np.dot(lambT, u2)[5], 0, 0]) / np.dot(lambT, u2)[5]))

print("norm")
print(np.linalg.norm(u2))

import pandas as pd

pd.DataFrame(lmid).to_csv(file_name, header=None, index=None, float_format='%.15f')

# In[13]:


# isocontour plot of DC potential from x (axial) direction
from electrode import System, GridElectrode

# strs = "DC1 DC2 DC3 DC4 DC5 DC6 DC7 DC8 DC9 DC10 DC11 DC12 DC13 DC14 DC15 DC16 DC17 DC18 DC19 DC20 DC21".split()
strs = "DC1 DC2 DC3 DC4 DC5 DC6 DC7 DC8 DC9 DC10 DC11 DC12 DC13 DC14 DC15 DC16 DC17 DC18 DC19 DC20 DC21".split()
Vx = np.zeros(())
# Vx = np.zeros((26,38))
i = 0
for inp in strs:
    arlo = np.array([])
    if inp in excl:
        print(inp)
    else:
        r = Result.from_vtk(prefix + suffix, inp)
        if inp == "RF":
            r = Result.from_vtk(prefix + suffix, "RF")
            p = r.pseudo_potential
        else:
            p = r.potential
        maxp = np.amax(p)
        #     x = grid.to_mgrid()[:,p.shape[0]//2]
        #     p = p[p.shape[0]//2]
        x = grid.to_mgrid()[:, :, p.shape[1] // 2]
        p = p[:, p.shape[1] // 2]
        if inp == "RF":
            print("trip")
            Vx = Vx + p * 0
        else:
            #         s[inp].dc = u2[i]*1
            Vx = Vx + (p * u2[i] * 1)
    #     else:
    #         Vx = Vx
    print(i)
    i = i + 1

# print("yz plane, %s potential"%ele)
fig, ax = plt.subplots()
ax.set_aspect("equal")
# yz plane should use x[1], x[2]. wwc
fig.set_size_inches(20, 10)
# ax.set_ylim(-1,2)
ax.grid(axis='both')
yticks = np.arange(-10, 10, 1)
ax.set_yticks(yticks)
xticks = np.arange(-10, 10, 1)
ax.set_xticks(xticks)
# ax.plot(3.7,1.8,marker = 'o')
ax.plot(xl, zl, marker='o', color='k')
ax.set_ylim(zl - Lz / 2, zl + Lz / 2)
ax.set_xlim(xl - Lx / 2, xl + Lx / 2)
print(Vx.min())
print(Vx.max())
print(np.shape(Vx))
yticks = np.arange(zl - Lz / 2, zl + Lz / 2, 0.02)
ax.set_yticks(yticks)
xticks = np.arange(xl - Lx / 2, xl + Lx / 2, 0.02)
ax.set_xticks(xticks)
# ax.contourf(x[1], x[2], Vx, levels=np.linspace(-10,10,20), cmap=plt.cm.RdYlGn)
# 2e-2
ax.contourf(x[0], x[2], Vx, levels=np.linspace(Vx.min(), Vx.max(), 100), cmap=plt.cm.RdYlGn)  # 2e-2

# In[198]:


# isocontour plot of DC potential from x (axial) direction
from electrode import System, GridElectrode

strs = "DC1 DC2 DC3 DC4 DC5 DC6 DC7 DC8 DC9 DC10 DC11 DC12 DC13 DC14 DC15 DC16 DC17 DC18 DC19 DC20 DC21 RF".split()
pseed = r.potential
Vx = np.zeros(())
i = 0
for inp in strs:
    if inp in excl:
        print("trigger")
    else:
        arlo = np.array([])
        r = Result.from_vtk(prefix + suffix, inp)
        if inp == "RF":
            p = r.pseudo_potential
        # method = 'Newton-CG'
        else:
            p = r.potential
        maxp = np.amax(p)
        x = grid.to_mgrid()[:, p.shape[0] // 2]
        p = p[p.shape[0] // 2]
        if inp == "RF":
            print("trip")
            Vx = Vx + p * 0
        else:
            val = p * u2[i] * 1
            #         s[inp].dc = u2[i]
            Vx = Vx + val
    #     else:
    #         Vx = Vx
    i = i + 1

# print("yz plane, %s potential"%ele)
fig, ax = plt.subplots()
ax.set_aspect("equal")
# yz plane should use x[1], x[2]. wwc
fig.set_size_inches(20, 10)
# ax.set_ylim(-1,2)
ax.grid(axis='both')
# ax.plot(3.7,1.8,marker = 'o')
ax.plot(yl, zl, marker='o', color='k')
ax.set_xlim(yl - Ly / 2, yl + Ly / 2)
ax.set_ylim(zl - Lz / 2, zl + Lz / 2)
xticks = np.arange(yl - Ly / 2, yl + Ly / 2, 0.02)
ax.set_xticks(xticks)
yticks = np.arange(zl - Lz / 2, zl + Lz / 2, 0.02)
ax.set_yticks(yticks)
print(Vx.min())
print(Vx.max())
print(np.shape(Vx))
# ax.contourf(x[1], x[2], Vx, levels=np.linspace(-10,10,20), cmap=plt.cm.RdYlGn)    # 2e-2
v = np.linspace(Vx.min(), Vx.max(), 80)
ax.contourf(x[1], x[2], Vx, levels=v, cmap=plt.cm.RdYlGn)  # 2e-2

# In[ ]:


# In[162]:


import scipy.constants as ct

l = 72e-6  # length scale
u = 100.  # peak rf voltage
o = 36e6 * 2 * np.pi  # rf frequency
m = 40 * ct.atomic_mass  # ion mass
q = 1 * ct.elementary_charge  # ion charge
# coefficient of pseudo-potential. See blakestad2010 Eq.(5.2). Run rf_scale() before other calculations.

rf_scale = s.rf_scale(m, q, l, o)
s["RF"].rf = u  # peak rf voltage
# method = 'Newton-CG'

# x0 = s.minimum((4.34507963, -0.04303287,  0.99403176))
x0 = np.array([3.89080547, -0.05081175, 1.05714154])
print(s.electrical_potential(x0)[0])
# for _ in s.analyze_static(x0,m=m,l=l, o=o):
#     print(_)

n = 30
# xyz = np.mgrid[-.1:.1:1j*n, -.1:.1:1j*n, 1.12:2]
# xyz = np.mgrid[0:1, -.02:.02:1j*n, .5:1.5:1j*n]
xyz = grid.to_mgrid()
p = s.potential(xyz.reshape(3, -1).T, 0).reshape(xyz[0].shape)
v = np.linspace(-3, 100, 17)
fig, ax = plt.subplots()
ax.set_aspect("equal")
fig.set_size_inches(8, 10)
ax.contourf(xyz[1, 10, :, :], xyz[2, 10, :, :], p[10, :, :], v, cmap=plt.cm.RdYlGn)


