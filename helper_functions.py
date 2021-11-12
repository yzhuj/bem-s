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
    fig, ax = plt.subplots(subplot_kw=dict(aspect="equal"), figsize=(12, 6), dpi=200)
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
    job.solve_singularities(num_mom=5, num_lev=1)
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
    x,y,z = grid.to_xyz() # p.shape[0]/2 is in the middle of x.
    p = p[p.shape[0] // 2,:,:] # get a slice of yz plane at x = p.shape[0]/2.
    print("yz plane, RF pseudo")
    fig, ax = plt.subplots()
    fig.set_size_inches(4, 10)
    ax.set_aspect("equal")
    ax.grid(axis='both')
    yticks = np.arange(0.5, 1.5, 0.1)
    ax.set_yticks(yticks)
    xticks = np.arange(-1, 1, 0.1)
    ax.set_xticks(xticks)
    # ax.set_ylim(0.5, 1.5)
    # ax.set_xlim(0.5,1.5)
    ax.contour(y,z, np.transpose(p), levels=np.linspace(p.min(),(p.max()-p.min())*0.1+p.min(), 100), cmap=plt.cm.RdYlGn)
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
    ax.contour(yp,xp, p, levels=np.linspace(p.min(),p.max(), 20), cmap=plt.cm.RdYlGn)
    plt.show()

def find_saddle(V,X,Y,Z,dim, scale=1, Z0=None,min=False):
    """Returns the indices of the local extremum or saddle point of the scalar A as (Is,Js,Ks).
    V is a 3D matrix containing an electric potential and must solve Laplace's equation
    X,Y,Z are the vectors that define the grid in three directions
    Z0: Z coordinate for saddle finding in a 2D potential slice
    For dim==2, the values of A are linearly extrapolated from [Z0] and [Z0]+1
    to those corresponding to Z0 and Ks is such that z[Ks]<Z0, z[Ks+1]>=Z0."""

    if (dim==2 and Z0==None):
        return 'z0 needed for evaluation'
    if dim==3:
        if len(V.shape)!=3:
            return('Problem with find_saddle.m dimensionalities.')
        f=V/float(np.amax(V)) # Normalize field
        [Ex,Ey,Ez]=np.gradient(f,abs(X[1]-X[0])/scale,abs(Y[1]-Y[0])/scale,abs(Z[1]-Z[0])/scale)
        [Ex2,Ey2,Ez2]=np.gradient(f,abs(X[1]-X[0])/scale,abs(Y[1]-Y[0])/scale,abs(Z[1]-Z[0])/scale,edge_order=2)# grid spacing is automatically consistent thanks to BEM-solver
        E=np.sqrt(Ex**2+Ey**2+Ez**2) # magnitude of gradient (E field)
        if min== True:
            E=f
        m=E[0,0,0]
        origin=[0,0,0]
        for i in range(E.shape[0]):
            for j in range(E.shape[1]):
                for k in range(E.shape[2]):
                    if E[i,j,k]<m:
                        # if Ex2[i,j,k]+Ey2[i,j,k]+Ez2[i,j,k]<0:
                        m=E[i,j,k]
                        origin=[i,j,k]
        if origin[0]==(0 or V.shape[0]):
            print('find_saddle: Saddle out of bounds in  x (i) direction.\n')
            return origin
        if origin[0]==(0 or V.shape[1]):
            print('find_saddle: Saddle out of bounds in  y (j) direction.\n')
            return origin
        if origin[0]==(0 or V.shape[2]):
            print('find_saddle: Saddle out of bounds in  z (k) direction.\n')
            return origin
    #################################################################################################
    if dim==2: # Extrapolate to the values of A at z0.
        V2=V
        if len(V.shape)==3:
            Ks=0 # in case there is no saddle point
            for i in range(len(Z)):
                if Z[i-1]<Z0 and Z[i]>=Z0:
                    Ks=i-1
                    if Z0<1:
                        Ks+=1
            Vs=V.shape
            if Ks>=len(Z):
                return('The selected coordinate is at the end of range.')
            v1=V[:,:,Ks]
            v2=V[:,:,Ks+1]
            V2=v1+(v2-v1)*(Z0-Z[Ks])/(Z[Ks+1]-Z[Ks])
        V2s=V2.shape
        if len(V2s)!=2: # Old: What is this supposed to check? Matlab code: (size(size(A2),2) ~= 2)
            return('Problem with find_saddle.py dimensionalities. It is {}.'.format(V2s))
        f=V2/float(np.max(abs(V2)))
        [Ex,Ey]=np.gradient(f,abs(X[1]-X[0]),abs(Y[1]-Y[0]))
        E=np.sqrt(Ex**2+Ey**2)
        m=float(np.min(E))
        mr=E[0,0]
        Is,Js=1,1 # in case there is no saddle
        for i in range(E.shape[0]):
            for j in range(E.shape[1]):
                if E[i,j]<mr:
                    mr=E[i,j]
                    Is,Js=i,j
        origin=[Is,Js,Ks]
        if Is==1 or Is==V.shape[0]:
            print('find_saddle: Saddle out of bounds in  x (i) direction.\n')
            return origin
        if Js==1 or Js==V.shape[1]:
            print('find_saddle: Saddle out of bounds in  y (j) direction.\n')
            return origin
    return origin
