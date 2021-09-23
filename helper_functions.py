import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from bem.formats import stl
import logging, os

def load_file(Mesh,Electrodes,prefix,scale,use_stl=True):
    if not use_stl:
        # load electrode faces from loops
        ele = Electrodes.from_trap(open("%s.ele" % prefix), scale)
        # initial triangulation, area 20, quiet
        mesh = Mesh.from_electrodes(ele)
        mpl.rcParams['lines.linewidth'] = 0.2
        mesh.triangulate(opts="a0.01q25.")
    else:
        # load electrode faces from colored stl
        # s_nta is intermediate processed stl file.
        s_nta = stl.read_stl(open("%s.stl" % prefix, "rb"))
        mpl.rcParams['lines.linewidth'] = 0.2
        print("Import stl:", os.path.abspath("./" + prefix + ".stl"), "\n")
        print("Electrode colors (numbers):\n")
        mesh = Mesh.from_mesh(stl.stl_to_mesh(*s_nta, scale=scale / 1e-3, rename={0: "DC21"}))
    return mesh,s_nta


def plot_mesh(xl,yl,mesh,scale):
    # Plot triangle meshes.
    fig, ax = plt.subplots(subplot_kw=dict(aspect="equal"), figsize=(24, 12), dpi=200)
    ax.set_xlabel("x/l", fontsize=10)
    ax.set_ylabel("y/l", fontsize=10)
    ax.text(-1.5, 7, "l = %d um" % (scale / 1e-6), fontsize=12)
    ax.plot(xl, yl, marker='o', color='k')
    # ax.grid(axis = 'both')
    yticks = np.arange(-100, 100, 2)
    ax.set_yticks(yticks)
    xticks = np.arange(-100, 100, 2)
    ax.set_xticks(xticks)
    mesh.plot(ax)
    plt.show()

# Define calculation function.
def run_job(args):
    # job is Configuration instance.
    job, grid, prefix = args
    # refine twice adaptively with increasing number of triangles, min angle 25 deg.
    # job.adapt_mesh(triangles=1e2, opts="q25Q")
    # job.adapt_mesh(triangles=1e3, opts="q25Q")
    # solve for surface charges
    job.solve_singularities(num_mom=5, num_lev=3)
#     print("done")
    # get potentials and fields
    result = job.simulate(grid, field=job.name=="RF", num_lev=4)    # For "RF", field=True computes the field.
    result.to_vtk(prefix)
    print("finished job %s" % job.name)
    return job.collect_charges()

def plot_RF(Result,prefix,suffix,grid):
    result = Result.from_vtk(prefix + suffix, "RF")
    p = result.pseudo_potential
    maxp = np.amax(p)
    print("p max", maxp)
    x = grid.to_mgrid()[:, p.shape[0] // 2]  # p.shape[0]/2 is in the middle of x.
    p = p[p.shape[0] // 2]  # get a slice of yz plane at x = p.shape[0]/2.
    print("yz plane, RF pseudo")
    fig, ax = plt.subplots()
    fig.set_size_inches(4, 10)
    ax.set_aspect("equal")
    ax.grid(axis='both')
    yticks = np.arange(0.5, 1.5, 0.2)
    ax.set_yticks(yticks)
    xticks = np.arange(-1, 1, 0.2)
    ax.set_xticks(xticks)
    # ax.set_ylim(0.5, 1.5)
    # ax.set_xlim(0.5,1.5)
    ax.contourf(x[1], x[2], p, levels=np.linspace(p.min(),(p.max()-p.min())*0.05+p.min(), 50), cmap=plt.cm.RdYlGn)
    plt.show()

def plot_DC(Result,prefix,suffix,grid,strs,dir='x'):
    p = np.zeros(Result.from_vtk(prefix + suffix, strs[0]).potential.shape)
    for em in strs:
        ele = em
        result = Result.from_vtk(prefix + suffix, ele)
        pmid = result.potential
        maxp = np.amax(p)
        #     print("p max", maxp)
        #     print(np.shape(Vx))
        p = p+pmid
    x,y,z = grid.to_xyz()
    if dir== 'x':
        p = p[p.shape[0] // 2,:,:]
        xp = y
        yp = z
    elif dir== 'y':
        p = p[:,p.shape[1] // 2,:]
        xp = x
        yp = z
    else:
        p = p[:,:,p.shape[2] // 2]
        xp = x
        yp = y
    print("yz plane, %s potential" % ele)
    fig, ax = plt.subplots()
    ax.set_aspect("equal")
    # yz plane should use x[1], x[2]. wwc
    X, Y = np.meshgrid(xp, yp)
    fig.set_size_inches(20, 10)
    ax.contourf(yp,xp, p, levels=np.linspace(p.min(),p.max(), 20), cmap=plt.cm.RdYlGn)
    plt.show()
