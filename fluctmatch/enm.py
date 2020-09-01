from os import path, mkdir
from collections import OrderedDict
from shutil import copyfile
from subprocess import check_call
import numpy as np
import pandas as pd
import MDAnalysis
from fluctmatch.miscell import get_patch, check_dir_exist_and_make
from fluctmatch.charmm import Script
from fluctmatch.sequence import sequences

cmm_root = '/home/yizaochen/PycharmProjects/connect_macro_micro'
d_atomcgtype = {'O1P': 'P', 'P': 'P', 'O2P': 'P', 'O5\'': 'P',
                'C5\'': 'P', 'O3\'': 'P', 'C4\'': 'S', 'O4\'': 'S',
                'C1\'': 'S', 'C2\'': 'S', 'C3\'': 'S', 'O2\'': 'S',
                'N1': 'B', 'C2': 'B', 'N2': 'B', 'O2': 'B', 'N3': 'B',
                'C4': 'B', 'O4': 'B', 'N4': 'B', 'C5': 'B', 'C5M': 'B',
                'C6': 'B', 'N6': 'B', 'O6': 'B', 'N7': 'B', 'C8': 'B',
                'N9': 'B', 'C7': 'B'}


class ENMAgent:
    def __init__(self, root, host, type_na):
        self.host = host
        self.type_na = type_na
        self.host_folder = path.join(root, host)
        self.na_folder = path.join(self.host_folder, type_na)
        self.input_folder = path.join(self.na_folder, 'input')
        self.charmminp_folder = path.join(self.na_folder, 'charmm_inp')
        self.charmmdat_folder = path.join(self.na_folder, 'charmm_dat')
        self.mode_traj_folder = path.join(self.na_folder, 'mode_traj')
        self.ic_folder = path.join(self.na_folder, 'ic')
        self.mat_folder = path.join(self.na_folder, 'ic_fluct_mat')
        self.rtficstr_folder = path.join(self.na_folder, 'rtf_ic_str')
        self.datafolder = path.join(self.na_folder, 'data')
        self.backupfolder = path.join(self.datafolder, 'backup')
        self.initialize_folders()
        
        self.seq1 = sequences[self.host][self.type_na]['guide']
        self.seq2 = sequences[self.host][self.type_na]['target']

        self.enmprm = path.join(self.datafolder, 'na_enm.prm')
        self.mode0dcd = path.join(self.mode_traj_folder, 'mode.0.dcd')
        self.crd = path.join(self.input_folder,
                             '{0}.nohydrogen.avg.crd'.format(self.type_na))

        self.u = None
        self.map = None
        self.inverse_map = None
        self.residues_map = None
        self.atomid_map = None
        self.atomid_map_inverse = None
        self.atomname_map = None
        self.strandid_map = None
        self.resid_map = None
        self.mass_map = None

        self.ics = dict()
        self.avgs = dict()

    def initialize_folders(self):
        for folder in [self.host_folder, self.na_folder, self.input_folder,
                       self.charmminp_folder, self.charmmdat_folder,
                       self.mode_traj_folder, self.ic_folder, self.mat_folder,
                       self.rtficstr_folder, self.datafolder, self.backupfolder]:
            check_dir_exist_and_make(folder)

    def check_avg_crd(self):
        if not path.isfile(self.crd): 
            print(f'{self.crd} does not exist.')
            return False
        else:
            print(f'{self.crd} exists.')
            return True
            
    def copy_avg_crd(self, old_folder):
        old_avg_crd = path.join(old_folder, self.host, self.type_na, 'input', 'heavyatoms', f'{self.type_na}.nohydrogen.avg.crd')
        copyfile(old_avg_crd, self.crd)
        print(f'cp {old_avg_crd} {self.crd}')
    
    def check_mode0dcd(self):
        if not path.isfile(self.mode0dcd): 
            print(f'{self.mode0dcd} does not exist.')
            return False
        else:
            print(f'{self.mode0dcd} exists.')
            return True
        
    def copy_mode0dcd(self, old_folder):
        old_mode0dcd = path.join(old_folder, self.host, self.type_na, 'input', 'heavyatoms', 'bdna+bdna.nohydrogen.fitavg.dcd') 
        copyfile(old_mode0dcd, self.mode0dcd)
        print(f'cp {old_mode0dcd} {self.mode0dcd}')

    def set_mda_universe(self):
        self.u = MDAnalysis.Universe(self.crd, self.crd)

    def set_required_d(self):
        self.map, self.inverse_map, self.residues_map, self.atomid_map,\
        self.atomid_map_inverse, self.atomname_map, self.strandid_map,\
        self.resid_map, self.mass_map = self.__build_map()

    def __build_map(self):
        d1 = dict()  # key: selction, value: cgname
        d2 = dict()  # key: cgname,   value: selection
        d3 = dict()
        d4 = dict()  # key: cgname, value: atomid
        d5 = dict()  # key: atomid, value: cgname
        d6 = dict()  # key: cgname, value: atomname
        d7 = dict()  # key: cgname, value: strand_id
        d8 = dict()  # key: cgname, value: resid
        d9 = dict()  # key: cgname, value: mass
        atomid = 1
        segid1 = self.u.select_atoms("segid STRAND1")
        d3['STRAND1'] = dict()
        for i, atom in enumerate(segid1):
            cgname = 'A{0}'.format(i+1)
            selection = get_selection(atom)
            d1[selection] = cgname
            d2[cgname] = selection
            if atom.resid not in d3['STRAND1']:
                d3['STRAND1'][atom.resid] = list()
            d3['STRAND1'][atom.resid].append(cgname)
            d4[cgname] = atomid
            d5[atomid] = cgname
            d6[cgname] = atom.name
            d7[cgname] = 'STRAND1'
            d8[cgname] = atom.resid
            d9[cgname] = atom.mass
            atomid += 1
        segid2 = self.u.select_atoms("segid STRAND2")
        d3['STRAND2'] = dict()
        for i, atom in enumerate(segid2):
            cgname = 'B{0}'.format(i+1)
            selection = get_selection(atom)
            d1[selection] = cgname
            d2[cgname] = selection
            if atom.resid not in d3['STRAND2']:
                d3['STRAND2'][atom.resid] = list()
            d3['STRAND2'][atom.resid].append(cgname)
            d4[cgname] = atomid
            d5[atomid] = cgname
            d6[cgname] = atom.name
            d7[cgname] = 'STRAND2'
            d8[cgname] = atom.resid
            d9[cgname] = atom.mass
            atomid += 1
        return d1, d2, d3, d4, d5, d6, d7, d8, d9

    def get_all_pairs(self, cutoff):
        result = list()
        for selection1, cgname1 in self.map.items():
            atoms = self.u.select_atoms("around {0} ({1})".format(cutoff, selection1))
            for atom in atoms:
                selection2 = get_selection(atom)
                cgname2 = self.map[selection2]
                if cgname1 == cgname2:
                    continue
                result.append(AtomPair(cgname1, cgname2, selection1, selection2))
        return list(OrderedDict.fromkeys(result))
        #  return list(set(result))

    def write_make_enm_crd_input(self):
        inp = path.join(self.charmminp_folder, 'make_enm_crd.inp')
        
        na = 'bdna'
        supplement1 = get_patch(self.seq1, 1)
        supplement2 = get_patch(self.seq2, 2)

        crd1 = path.join(self.mkcrd_folder, '{0}1.crd'.format(na))
        inp1 = Script(path.join(self.mkcrd_folder, '{0}1.inp'.format(na)))
        inp1.write_bomlev()
        inp1.initialize_rtf_prm(amber=amber)
        inp1.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        if supplement1 is not None:
            inp1.write_supplement(supplement1)
        inp1.gen_angle_dihedral()
        inp1.read_pdb(path.join(self.mkcrd_folder, '{0}1.1.pdb'.format(na)))
        inp1.write_crd(crd1)
        inp1.end()

    def make_enm_crd(self):
        charmm = "/home/yizaochen/opt/charmm/exec/gnu/charmm"
        inp = path.join(self.charmminp_folder, 'make_enm_crd.inp')
        dat = path.join(self.charmmdat_folder, 'make_enm_crd.dat')
        self.__exec_charmm(charmm, inp, dat)
            
    def __exec_charmm(self, charmm, f_input, f_output):
        print("charmm< {0} > {1}".format(f_input, f_output))
        check_call(charmm, stdin=open(f_input, 'r'), stdout=open(f_output, 'w+'), shell=True)


    def read_ic_fluct_matrix(self, modeid):
        f_in = path.join(self.mat_folder, 'mode.{0}.npy'.format(modeid))
        return np.load(f_in)

    def write_ic_fluct_matrix(self, modeid):
        n_atoms = len(self.u.atoms)
        mat = np.zeros((n_atoms, n_atoms))
        data = self.get_ic_fluct(modeid)
        for pair in data:
            i = self.atomid_map[pair.name1] - 1
            j = self.atomid_map[pair.name2] - 1
            value = pair.value / 2
            mat[i, j] = value
            mat[j, i] = value
        f_out = path.join(self.mat_folder, 'mode.{0}.npy'.format(modeid))
        np.save(f_out, mat)
        return mat

    def get_ic_fluct(self, modeid):
        try:
            return self.ics[modeid]
        except KeyError:
            self.read_mode_fluct(modeid)
            return self.ics[modeid]

    def get_ic_avg(self, modeid):
        try:
            return self.avgs[modeid]
        except KeyError:
            self.read_mode_avg(modeid)
            return self.avgs[modeid]

    def get_fluct_mean_std(self, modeid):
        data = self.get_ic_fluct(modeid)
        temp = [pair.value for pair in data]
        temp = np.array(temp)
        return temp.mean(), temp.std()

    def get_flucts_larger_than_assign(self, modeid, value):
        data = self.get_ic_fluct(modeid)
        temp = [pair for pair in data if pair.value >= value]
        return temp

    def get_statistics(self, modeid):
        data = self.get_ic_fluct(modeid)
        d = {'fluct': [pair.value for pair in data]}
        df = pd.DataFrame(d)
        print(df.describe())

    def get_mapped_df(self, modeid, threshold=0):
        data = self.get_flucts_larger_than_assign(modeid, threshold)
        result = {'Site1': list(), 'Site2': list(), 'fluctuation': list()}
        for pair in data:
            result['Site1'].append(self.inverse_map[pair.name1])
            result['Site2'].append(self.inverse_map[pair.name2])
            result['fluctuation'].append(pair.value)
        df = pd.DataFrame(result)
        return df

    def get_df(self, modeid, threshold=0):
        data = self.get_flucts_larger_than_assign(modeid, threshold)
        result = {'Site1': list(), 'Site2': list(), 'fluctuation': list()}
        for pair in data:
            result['Site1'].append(pair.name1)
            result['Site2'].append(pair.name2)
            result['fluctuation'].append(pair.value)
        df = pd.DataFrame(result)
        return df

    def get_residues_decomposition(self, modeid):
        mat = self.read_ic_fluct_matrix(modeid)
        d = dict()
        summation = mat.sum()
        for strand in ['STRAND1', 'STRAND2']:
            d[strand] = dict()
            for resid in range(1, 11):
                sites = self.residues_map[strand][resid]
                atomids = get_atomids_from_sites(sites, self.atomid_map)
                min_id = min(atomids) - 1
                max_id = max(atomids) - 1
                sub_mat = mat[min_id:max_id+1, :]
                result = sub_mat.sum()
                d[strand][resid] = (result / summation) * 100
        return d

    def get_atoms_decomposition(self, modeid, strand, resid):
        mat = self.read_ic_fluct_matrix(modeid)
        sites = self.residues_map[strand][resid]
        atomids = get_atomids_from_sites(sites, self.atomid_map)

        min_id = min(atomids) - 1
        max_id = max(atomids) - 1
        sub_mat = mat[min_id:max_id + 1, :]
        summation = sub_mat.sum()

        d_result = dict()
        for site in sites:
            atomid = self.atomid_map[site]
            temp_mat = mat[atomid - 1, :]
            result = temp_mat.sum()
            d_result[site] = (result / summation) * 100

        d_result_1 = dict()
        for key, value in d_result.items():
            newkey = self.atomname_map[key]
            d_result_1[newkey] = value
        return d_result_1

    def get_atoms_decomposition_divide_n_pair(self, modeid, strand, resid):
        mat = self.read_ic_fluct_matrix(modeid)
        sites = self.residues_map[strand][resid]

        d_result = dict()
        l_result = list()
        for site in sites:
            atomid = self.atomid_map[site]
            temp_mat = mat[atomid - 1, :]
            n_pairs = len(np.nonzero(temp_mat)[0])
            result = temp_mat.sum() / n_pairs
            d_result[site] = result
            l_result.append(result)
        min_element = min(l_result)
        max_element = max(l_result)
        interval = max_element - min_element

        d_result_1 = dict()
        for key, value in d_result.items():
            newkey = self.atomname_map[key]
            d_result_1[newkey] = (value - min_element) / interval
        return d_result_1

    def read_mode_fluct(self, modeid):
        result = list()
        f_in = path.join(self.ic_folder, 'mode.{0}.modified.ic'.format(modeid))
        data = np.genfromtxt(f_in, dtype=str, skip_header=5)
        for subdata in data:
            result.append(FluctPair(subdata[2], subdata[4], float(subdata[9])))
            result.append(FluctPair(subdata[6], subdata[8], float(subdata[13])))
        self.ics[modeid] = list(set(result))

    def read_mode_avg(self, modeid):
        result = list()
        f_in = path.join(self.ic_folder, 'mode.{0}.avg.modified.ic'.format(modeid))
        data = np.genfromtxt(f_in, dtype=str, skip_header=5)
        for subdata in data:
            result.append(FluctPair(subdata[2], subdata[4], float(subdata[9])))
            result.append(FluctPair(subdata[6], subdata[8], float(subdata[13])))
        self.avgs[modeid] = list(set(result))

    def read_forceconstants(self):
        f_in = path.join(self.datafolder, 'na_enm.prm')
        data = np.genfromtxt(f_in, skip_header=4, skip_footer=2, dtype=str)
        d_result = {'Site1': list(), 'Site2': list(), 'k': list(), 'b': list(),
                    'Name1': list(), 'Name2': list()}
        for subdata in data:
            d_result['Site1'].append(subdata[0])
            d_result['Site2'].append(subdata[1])
            d_result['k'].append(float(subdata[2]))
            d_result['b'].append(float(subdata[3]))
            d_result['Name1'].append(self.inverse_map[subdata[0]])
            d_result['Name2'].append(self.inverse_map[subdata[1]])
        df = pd.DataFrame(d_result)
        df = df[['Name1', 'Name2', 'k', 'b']]
        return df

    def read_forceconstants_site_name(self):
        f_in = path.join(self.datafolder, 'na_enm.prm')
        data = np.genfromtxt(f_in, skip_header=4, skip_footer=2, dtype=str)
        d_result = {'Site1': list(), 'Site2': list(), 'k': list(), 'b': list(),
                    'Name1': list(), 'Name2': list()}
        for subdata in data:
            d_result['Site1'].append(subdata[0])
            d_result['Site2'].append(subdata[1])
            d_result['k'].append(float(subdata[2]))
            d_result['b'].append(float(subdata[3]))
            d_result['Name1'].append(self.inverse_map[subdata[0]])
            d_result['Name2'].append(self.inverse_map[subdata[1]])
        df = pd.DataFrame(d_result)
        df = df[['Site1', 'Site2', 'k', 'b']]
        return df

    def read_forceconstants_accord_type(self):
        f_in = path.join(self.datafolder, 'na_enm.prm')
        data = np.genfromtxt(f_in, skip_header=4, skip_footer=1, dtype=str)
        columns = ['Strand_i', 'Resid_i', 'Name_i', 'Type_i',
                   'Strand_j', 'Resid_j', 'Name_j', 'Type_j',
                   'k']
        d_result = dict()
        for name in columns:
            d_result[name] = list()
        for subdata in data:
            d_result['Strand_i'].append(self.strandid_map[subdata[0]])
            d_result['Strand_j'].append(self.strandid_map[subdata[1]])
            d_result['Resid_i'].append(self.resid_map[subdata[0]])
            d_result['Resid_j'].append(self.resid_map[subdata[1]])
            d_result['Name_i'].append(self.atomname_map[subdata[0]])
            d_result['Name_j'].append(self.atomname_map[subdata[1]])
            d_result['Type_i'].append(d_atomcgtype[self.atomname_map[subdata[0]]])
            d_result['Type_j'].append(d_atomcgtype[self.atomname_map[subdata[1]]])
            d_result['k'].append(float(subdata[2]))
        df = pd.DataFrame(d_result)
        df = df[columns]
        return df

    def read_k_dx_mode(self, modeid):
        f_in = path.join(self.datafolder, 'na_enm.prm')
        data = np.genfromtxt(f_in, skip_header=4, skip_footer=1, dtype=str)
        columns = ['Strand_i', 'Resid_i', 'Name_i', 'Type_i',
                   'Strand_j', 'Resid_j', 'Name_j', 'Type_j',
                   'k', 'dx']
        d_result = dict()
        for name in columns:
            d_result[name] = list()
        for subdata in data:
            d_result['Strand_i'].append(self.strandid_map[subdata[0]])
            d_result['Strand_j'].append(self.strandid_map[subdata[1]])
            d_result['Resid_i'].append(self.resid_map[subdata[0]])
            d_result['Resid_j'].append(self.resid_map[subdata[1]])
            d_result['Name_i'].append(self.atomname_map[subdata[0]])
            d_result['Name_j'].append(self.atomname_map[subdata[1]])
            d_result['Type_i'].append(d_atomcgtype[self.atomname_map[subdata[0]]])
            d_result['Type_j'].append(d_atomcgtype[self.atomname_map[subdata[1]]])
            d_result['k'].append(float(subdata[2]))
        f_in = path.join(self.ic_folder, 'mode.{0}.modified.ic'.format(modeid))
        data = np.genfromtxt(f_in, dtype=str, skip_header=5)
        for subdata in data:
            d_result['dx'].append(float(subdata[9]))
        df = pd.DataFrame(d_result)
        df = df[columns]
        return df

    def get_avg_structure(self):
        d = read_structure(self.crd)
        return d


class AtomPair:
    def __init__(self, name1, name2, selection1, selection2, k=10., b=5.):
        self.name1 = name1
        self.name2 = name2
        self.selection1 = selection1
        self.selection2 = selection2
        self.atom1 = None  # MDAnalysis Atom Obejct
        self.atom2 = None  # MDAnalysis Atom Obejct
        self.distances = list()
        self.k = k
        self.b = b
        self.pairname = "{0}-{1}".format(self.name1, self.name2)

    def set_atom1(self, atom):
        self.atom1 = atom

    def set_atom2(self, atom):
        self.atom2 = atom

    def get_distance(self):
        return np.linalg.norm(self.atom1.positions[0] - self.atom2.positions[0])

    def append_distances(self):
        self.distances.append(np.linalg.norm(self.atom1.positions[0] - self.atom2.positions[0]))

    def __repr__(self):
        return "{0}-{1}".format(self.name1, self.name2)

    def __hash__(self):
        return hash(self.name1) + hash(self.name2)

    def __eq__(self, other):
        if self.name1 == other.name1 and self.name2 == other.name2:
            return True
        elif self.name1 == other.name2 and self.name2 == other.name1:
            return True
        else:
            return False


class FluctPair(AtomPair):
    def __init__(self, name1, name2, fluct_value):
        self.name1 = name1
        self.name2 = name2
        self.value = fluct_value


def check_dir_exist_and_make(file_path):
    if path.exists(file_path):
        print("{0} exists".format(file_path))
    else:
        print("mkdir {0}".format(file_path))
        mkdir(file_path)


def get_selection(atom):
    return 'segid {0} and resid {1} and name {2}'.format(atom.segid, atom.resid, atom.name)


def df_2_fluctpairs(df):
    l_store = list()
    for index, row in df.iterrows():
        pair = FluctPair(row['Site1'], row['Site2'], row['fluctuation'])
        l_store.append(pair)
    return l_store


def get_atomids_from_sites(sites, d_map):
    l_store = [d_map[site] for site in sites]
    return l_store


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


def read_hessian(f_in, natoms):
    # Read File
    f = open(f_in, 'r')
    lines = f.readlines()
    f.close()

    # Useful Variables
    n_3 = 3 * natoms  # 3N
    nij = int(n_3 * (n_3 + 1) / 2)  # The number of elements of upper triangle matrix
    first_line = 2 + natoms
    last_line = first_line + nij

    # Using genfromtxt to quickly make array
    temp = np.genfromtxt(lines[first_line:last_line])

    # Make Symmetry Matrix
    hessian = np.zeros((n_3, n_3))
    i_upper = np.triu_indices(n_3)  # indices of upper triangle contains diagonal
    i_lower = np.tril_indices(n_3, -1)  # indices of lower triangle without diagonal
    hessian[i_upper] = temp
    hessian[i_lower] = hessian.T[i_lower]
    return hessian


def create_upper_matrix(values, size):
    upper = np.zeros((size, size))
    upper[np.triu_indices(size, 0)] = values
    return upper

