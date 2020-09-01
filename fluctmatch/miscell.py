from os import path, mkdir
import pickle
import numpy as np

def get_patch(seq, strand_id):
    patchlist = list()
    terminate = len(seq) - 1
    for i, nt in enumerate(seq):
        if nt == 'T':
            continue
        if i == 0:
            temp = f'patch DO{nt}5 strand{strand_id} {i+1}\n'
        elif i == terminate:
            temp = f'patch DO{nt}3 strand{strand_id} {i+1}\n'
        else:
            temp = f'patch DO{nt} strand{strand_id} {i+1}\n'
        patchlist.append(temp)
    return ''.join(patchlist)


def get_patch_for_dna():
    patch1 = list()
    patch2 = list()


def get_antistrand_resid(guide_resid, n_bp=10):
    # This Function is ad-hoc, it should be changed for different systems
    return n_bp - guide_resid + 1


def check_dir_exist_and_make(file_path):
    if path.exists(file_path):
        print("{0} exists".format(file_path))
    else:
        print("mkdir {0}".format(file_path))
        mkdir(file_path)


def check_file_exist(file_path):
    exists = path.isfile(file_path)
    if exists:
        print("{0} already exist".format(file_path))
        return True
    else:
        print("{0} not exist!!!!!!".format(file_path))
        return False


def read_structure(f_in):
    d = dict()
    data = np.genfromtxt(f_in, skip_header=4)
    for subdata in data:
        key = int(subdata[0])
        d[key] = dict()
        d[key]['x'] = subdata[4]
        d[key]['y'] = subdata[5]
        d[key]['z'] = subdata[6]
    return d


def get_vector(f_in, f_format='xyz'):
    if f_format == 'crd':
        d = read_structure(f_in)
    elif f_format == 'xyz':
        d = read_structure_from_xyz(f_in)
    n_atoms = len(d)
    result = list()
    for i in range(1, n_atoms+1):
        for j in ['x', 'y', 'z']:
            result.append(d[i][j])
    return np.array(result)


def get_vector_by_d(d):
    n_atoms = len(d)
    result = list()
    for i in range(1, n_atoms+1):
        for j in ['x', 'y', 'z']:
            result.append(d[i][j])
    return np.array(result)


def get_d_by_vector(vector):
    d = dict()
    atomid = 0
    for idx, data in enumerate(vector):
        if idx % 3 == 0:
            atomid += 1
            d[atomid] = dict()
            d[atomid]['x'] = data
        elif idx % 3 == 1:
            d[atomid]['y'] = data
        else:
            d[atomid]['z'] = data
    return d


def read_structure_from_xyz(f_in):
    d = dict()
    data = np.genfromtxt(f_in, skip_header=1)
    for subdata in data:
        key = int(subdata[0])
        d[key] = dict()
        d[key]['x'] = subdata[1]
        d[key]['y'] = subdata[2]
        d[key]['z'] = subdata[3]
    return d


def get_distance_betw2atoms_from_d(d, atomid1, atomid2):
    v1 = [d[atomid1]['x']-d[atomid2]['x'], d[atomid1]['y']-d[atomid2]['y'], d[atomid1]['z']-d[atomid2]['z']]
    return np.linalg.norm(v1)


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            # >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            # >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            # >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def get_modulus_angle_between_two_vectors(v1, v2):
    v1_modulus = np.linalg.norm(v1)
    v2_modulus = np.linalg.norm(v2)
    v1_u = v1 / v1_modulus
    v2_u = v2 / v2_modulus
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    return v1_modulus, v2_modulus, angle


def rotation_mat_between_2vect(vector_orig, vector_fin):
    """Calculate the rotation matrix required to rotate from one vector to another.
    For the rotation of one vector to another, there are an infinit series of rotation matrices
    possible.  Due to axially symmetry, the rotation axis can be any vector lying in the symmetry
    plane between the two vectors.  Hence the axis-angle convention will be used to construct the
    matrix with the rotation axis defined as the cross product of the two vectors.  The rotation
    angle is the arccosine of the dot product of the two unit vectors.
    Given a unit vector parallel to the rotation axis, w = [x, y, z] and the rotation angle a,
    the rotation matrix R is::
              |  1 + (1-cos(a))*(x*x-1)   -z*sin(a)+(1-cos(a))*x*y   y*sin(a)+(1-cos(a))*x*z |
        R  =  |  z*sin(a)+(1-cos(a))*x*y   1 + (1-cos(a))*(y*y-1)   -x*sin(a)+(1-cos(a))*y*z |
              | -y*sin(a)+(1-cos(a))*x*z   x*sin(a)+(1-cos(a))*y*z   1 + (1-cos(a))*(z*z-1)  |
    @param vector_orig: The unrotated vector defined in the reference frame.
    @type vector_orig:  numpy array, len 3
    @param vector_fin:  The rotated vector defined in the reference frame.
    @type vector_fin:   numpy array, len 3
    """

    # Convert the vectors to unit vectors.
    vector_orig = vector_orig / np.linalg.norm(vector_orig)
    vector_fin = vector_fin / np.linalg.norm(vector_fin)

    # The rotation axis (normalised).
    axis = np.cross(vector_orig, vector_fin)
    axis_len = np.linalg.norm(axis)
    if axis_len != 0.0:
        axis = axis / axis_len

    # Alias the axis coordinates.
    x = axis[0]
    y = axis[1]
    z = axis[2]

    # The rotation angle.
    angle = np.arccos(np.dot(vector_orig, vector_fin))

    # Trig functions (only need to do this maths once!).
    ca = np.cos(angle)
    sa = np.sin(angle)

    # Initialize Rotation Matrix
    rot_mat = np.zeros((3, 3))

    # Calculate the rotation matrix elements.
    rot_mat[0, 0] = 1.0 + (1.0 - ca)*(x**2 - 1.0)
    rot_mat[0, 1] = -z*sa + (1.0 - ca)*x*y
    rot_mat[0, 2] = y*sa + (1.0 - ca)*x*z
    rot_mat[1, 0] = z*sa+(1.0 - ca)*x*y
    rot_mat[1, 1] = 1.0 + (1.0 - ca)*(y**2 - 1.0)
    rot_mat[1, 2] = -x*sa+(1.0 - ca)*y*z
    rot_mat[2, 0] = -y*sa+(1.0 - ca)*x*z
    rot_mat[2, 1] = x*sa+(1.0 - ca)*y*z
    rot_mat[2, 2] = 1.0 + (1.0 - ca)*(z**2 - 1.0)
    return rot_mat


def get_rotation_matrix(angle, axis='x'):
    ca = np.cos(angle)
    sa = np.sin(angle)
    if axis == 'x':
        rot_mat = [[1, 0, 0],
                   [0, ca, -sa],
                   [0, sa, ca]]
    elif axis == 'y':
        rot_mat = [[ca, 0, sa],
                   [0, 1, 0],
                   [-sa, 0, ca]]
    else:
        rot_mat = [[ca, -sa, 0],
                   [sa, ca, 0],
                   [0, 0, 1]]
    return rot_mat


def get_bar_x_y(lastmode, mat):
    y = list()
    for mode in range(1, lastmode+1):
        y.append(mat[mode-1, mode-1])
    x = np.arange(len(y))
    return x, y


def get_position_array(structure, atomid):
    return np.array([structure[atomid]['x'], structure[atomid]['y'], structure[atomid]['z']])


def save_dict_to_pkl(d, f_out):
    f = open(f_out, 'wb')
    pickle.dump(d, f)
    f.close()


def load_pkl_from_dict(f_in):
    f = open(f_in, 'rb')
    data = pickle.load(f)
    f.close()
    return data


