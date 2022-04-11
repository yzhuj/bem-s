import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from bem.formats import stl
import os
from scipy.signal import argrelextrema
from bem import Result
import pickle

def load_file(Mesh,Electrodes,prefix,scale,use_stl=True):
    if not use_stl:
        # load electrode faces from loops
        ele = Electrodes.from_trap(open("%s.ele" % prefix), scale)
        # initial triangulation, area 20, quiet
        mesh = Mesh.from_electrodes(ele)
        mesh.triangulate(opts="a0.01q25.")
    else:
        # load electrode faces from colored stl
        # s_nta is intermediate processed stl file.
        s_nta = stl.read_stl(open("trapstl/%s.stl" % prefix, "rb"))
        mesh = Mesh.from_mesh(stl.stl_to_mesh(*s_nta, scale=scale / 1e-3, rename={0: "DC21"}))
    return mesh,s_nta

#Create custom subplot w/ dimensions that you want, add the trapping point,
#then use mesh object's 'plot' function to add the mesh to it.
def plot_mesh(xl,yl,mesh,scale):
    # Plot triangle meshes.
    fig, ax = plt.subplots(subplot_kw=dict(aspect="equal"), figsize=(12, 6), dpi=400)
    ax.set_xlabel("x/l", fontsize=10)
    ax.set_ylabel("y/l", fontsize=10)
    ax.text(0, 0, "l = %d um" % (scale / 1e-6), fontsize=12)
    ax.plot(xl, yl, marker='.', color='k')
    # ax.grid(axis = 'both')
    yticks = np.arange(-1, 1, 0.1)
    ax.set_yticks(yticks)
    xticks = np.arange(-1, 1, 0.1)
    ax.set_xticks(xticks)
    mesh.plot(ax)
    plt.show()

# Trap simulations.
def run_job(args):
    # job is Configuration instance.
    job, grid, prefix = args

    # refine twice adaptively with increasing number of triangles, min angle 25 deg.
    # job.adapt_mesh(triangles=4e2, opts="q25Q")
    # job.adapt_mesh(triangles=1e4, opts="q25Q")
    # solve for surface charges
    job.solve_singularities(num_mom=4, num_lev=3)
#     print("done")
    # get potentials and fields
    # For "RF", field=True computes the field.
    result = job.simulate(grid, field=job.name=="RF", num_lev=2)
    result.to_vtk(prefix)
    print("finished job %s" % job.name)
    return job.collect_charges()

def write_pickle(fin,fout,grid):
    #grid is the field grid pts that give the locations of each simulated potential point
    #fin is the filename of the of the input vtk sim file
    #fout is the filename of the pickle you want to save to
    x, y, z = grid.to_xyz()
    nx = len(x)
    ny = len(y)
    nz = len(z)
    ntotal = nx * ny * nz

    trap = {'X': x,
            'Y': y,
            'Z': z}
    i = 0
    strs = "DC1 DC2 DC3 DC4 DC5 DC6 DC7 DC8 DC9 DC10 DC11 DC12 DC13 DC14 DC15 DC16 DC17 DC18 DC19 DC20 DC21".split()
    result0 = Result.from_vtk(fin, 'DC1')
    p0 = result0.potential
    for ele in strs:
        # if ele not in excl:
        result = Result.from_vtk(fin, ele)
        p = result.potential
        p = np.swapaxes(p, 0, 2)
        p = np.swapaxes(p, 0, 1)
        trap[ele] = {'potential': p}
        trap[ele]['position'] = [0, i]
        # else:
        #     trap[ele] = {'potential': np.zeros(np.shape(p0))}
        #     trap[ele]['position'] = [0, i]
        # i = i + 1

    electrode_list = strs

    f = open('./'+fout+'.pkl', 'wb')
    trap1 = {'X': trap['Y'],
             'Y': trap['Z'],
             'Z': trap['X'],
             'electrodes': {}}
    for electrode in electrode_list:
        trap1['electrodes'][electrode] = trap[electrode]
    pickle.dump(trap1, f, -1)
    f.close()