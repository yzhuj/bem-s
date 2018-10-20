from collections import OrderedDict
import struct

import numpy as np


def read_stl(fil, ignore_attr_length=True):     # wwc stl binary format
    """reads binary STL in file "fil" with Inventor style
    (color hidden in attribute length), returns array of normals, array
    of vertices and list or attributes

    Inventor >= 2013 supports colors (bottom right in the export options
    dialog)
    
    The atribute short in binary STL is used by some tools to hide
    color information. The format is 15 bit RGB or BGR, MSB is validity
    (default color if MSB is 0) COLOR=\\xbf\\xbf\\xbf\\xff in the header
    is the overall color default color.
    """
    h = fil.read(80)        # wwc A ASCII header of length 80 bytes.
    assert h.startswith("STLB")     # wwc see if it's a binary stl or raise an AssertionError.
    n, = struct.unpack("< I", fil.read(4))      # wwc n may be number of faces in the file. 4 bytes unsigned int. "< I" google python struct.unpack and see https://docs.python.org/3.6/library/struct.html#format-strings
    if ignore_attr_length:
        dtype = np.dtype([
            ("normal", "<f4", (3,)),
            ("triangle", "<f4", (3, 3)),
            ("attribute", "<u2")])
        data = np.fromfile(fil, dtype, n)
        normals, triangles = data["normal"], data["triangle"]
        attributes = data["attribute"]
    else: # actual STL binary format with variable length attribute bytes
        vectors, attributes = [], []
        while True:
            try:
                vectors.append(struct.unpack("< 12f", fil.read(48)))
                na, = struct.unpack("< H", fil.read(2))
                attributes.append(fil.read(na))
            except struct.error:
                break
        assert len(vectors) == n, (len(vectors), n)
        vectors = np.array(vectors)
        normals = vectors[:, :3]
        triangles = vectors[:, 3:].reshape(-1, 3, 3)
    return normals, triangles, attributes


def check_normals(normals, triangles):
    """verifies that given vertices are right-handed around normals"""
    a, b, c = triangles.transpose((1, 0, 2))
    n = np.cross(b-a, c-a)
    n /= np.sqrt(np.sum(np.square(n), axis=1))[:, None]
    assert np.allclose(n, normals, rtol=1e-3, atol=1e-10)

# wwc function
def partition_normals(normals,triangles,numbers=[],TOL=1e-6):
    """Partition points into different planes according to the normals in 
    one electrode, for 3D/2D meshing (Shewchuk's triangle C code is 2D meshing).
    Return points_numbers which has following format

    [(plane1_points,serial_numbers1),(plane2_points,serial_numbers2),...]

    plane_points -- array([ [x1,y1,z1],[x2,y2,z2],... [xn,yn,zn] ])
    serial_numbers -- array([ [0,1,2],[3,4,5],... [n-3,n-2,n-1] ])

    TOL: normal deviation tolerances
    """
    nm_unique = [normals[0]]      # normals[0] as first nm  wwc
    points = [[]]    # Each sublist [] represents a plane  wwc
    for nm, tr in zip(normals, triangles):
        add_plane = True
        # Search existing normals in nm_unique, in reversed order (points in the same plane are often ajacent.)  wwc
        for ith in range(len(nm_unique)-1,-1,-1):   # ith normal -- ith plane  wwc
            if np.linalg.norm(nm-nm_unique[ith]) < TOL:
                # points[ith] plane -- [array([x1,y1,z1]),array([x2,y2,z2]), ...]  wwc
                points[ith].extend(tr)
                add_plane = False
                break
        if add_plane:   # All normals don't coincide with the new nm, record it.  wwc
            nm_unique.append(nm)
            points.append([])
            points[-1].extend(tr)
    points_numbers = []
    for plane in points:
        plane_points = np.ascontiguousarray(plane,dtype=np.double)
        i = np.arange(0,plane_points.shape[0],3)    # shape[0] is the number of points in the plane.  wwc
        serial_numbers = np.c_[i, i+1, i+2].astype(np.intc)
        points_numbers.append( (plane_points,serial_numbers) )
    return points_numbers, np.array(nm_unique,dtype=np.double)

# wwc function. Copied from cpy.py. This function can split normals merely by 3 points coordinates of triangles, 
# though we alreadly have normals imported from stl file. So use above partition_normals() is OK.  wwc
def split_by_normal(stlmesh):
    """split curved faces into one face per triangle (aka split by
    normal, planarize). in place"""
    for name, faces in stlmesh.iteritems():    # An iterator  wwc
        new_faces = []
        # faces = [(points,triangles)], so points, triangles = (points, triangles)  wwc
        for points, triangles in faces:
            x = points[triangles, :]    # triangles are the serial number of points.  wwc
            # For all triangle points x0,x1,x2, get the cross product of vector x1-x0, x2-x0, the normal vectors.  wwc
            normals = np.cross(x[:, 1]-x[:, 0], x[:, 2]-x[:, 0])
            normals /= np.sqrt(np.sum(np.square(normals), axis=1))[:, None]
            if np.allclose(normals, normals[0][None, :]):
                new_faces.append((points, triangles))
            else:
                for triangle in triangles:
                    new_faces.append((points[triangle, :],
                        np.arange(3, dtype=np.intc).reshape((1, 3))))
        stlmesh[name] = new_faces
    return stlmesh


def stl_to_mesh(normals, triangles, attributes, scale=None, rename=None):
    """generates a {name: [(points, triangles)]} mesh from the
    stl arrays. For further treatment, use e.g:

    >>> s = read_stl(open(filename, "rb"))
    >>> check_normals(*s[:2])
    >>> r = stl_to_mesh(*s)
    >>> del r["stl_0"] # delete the color-less faces
    >>> print r.keys()
    ['stl_18761', 'stl_12943']
    >>> r = split_by_normal(r)
    >>> m = Mesh.from_mesh(r)
    >>> m.triangulate("qQ")
    >>> m.to_vtk("test_stl_import")
    
    rename can be a renaming dictionary mapping stl color numbers to
    electrode names. if None the color name is "stl_%i" %
    color_numer, else if the color is not found in rename,
    the electrode is dropped.
    """
    d = OrderedDict()    # d is an internediate dict which will pass the items to o.  wwc
    for a, nm, tr in zip(attributes, normals, triangles):
        d.setdefault(a, [[],[]])[0].append(nm) 
        d[a][1].append(tr)
        # d[a] = [nms, trs], nms = [array([x,y,z]),...] trs = [ array([[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]]),...]  wwc
    o = OrderedDict()    # o stores the {name: [(points, triangles)]}  wwc
    for a, nm_tr in d.iteritems():    # a is stl color number.  wwc
        nms, trs = np.array(nm_tr[0]), np.array(nm_tr[1])
        if scale:    # scale a dimensionless number. For example, scale = 40.  wwc
            # points /= scale
            trs /= scale    # nms are unit vectors.  wwc
        partition = False
        if rename:
            if a not in rename:
                # If you don't know the number of color a, this prints it for you to rename.  wwc
                print "dropping", a
                continue
            else:
                n = rename[a]    # If you have provided a rename for a, this renames it.
                partition = True
        else:
            n = "stl_%i" % a
        if partition:    # Designed for 3D multiplane partition.  wwc  
            o[n], planes = partition_normals(nms, trs)
            print "%i planes in electrode %s"%(len(planes),n)
            print "Normals vectors are:\n", planes
        else:    # Robert's original 2D version.  wwc
            i = np.arange(0, trs.shape[0]*3, 3)
            points = trs.reshape(-1, 3).astype(np.double)    # reshape puts all points together, not distinguish different triangles.  wwc
            triangles = np.c_[i, i+1, i+2].astype(np.intc)    # "triangles" here is actually serial number of points.  wwc
            # Note there is an additional list [ ] contains the tuple (points, triangles).  wwc
            # Robert might leave this as a future interface for 3D meshing.  wwc
            o[n] = [(np.ascontiguousarray(points),    # n = o.keys(), is "DC1","RF"..., the rename.  wwc
                np.ascontiguousarray(triangles))]     # triangles are the serial number of points.  wwc
        # o[n] = [(np.ascontiguousarray(points),
        #         np.ascontiguousarray(triangles))]
    return o


if __name__ == "__main__":
    import sys
    r = read_stl(open(sys.argv[1], "rb"))
    m = stl_to_mesh(*r)
    print m