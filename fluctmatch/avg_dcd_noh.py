from os import path
from shutil import copyfile
from subprocess import check_call
from charmm import Script
from miscell import check_dir_exist_and_make, get_patch
from sequence import sequences
import PDB

charmm = '/home/yizaochen/opt/charmm/exec/gnu/charmm'


class AvgcrddcdAgent:
    def __init__(self, host, type_na, rootfolder):
        self.host = host
        self.type_na = type_na
        self.host_folder = path.join(rootfolder, host)
        self.na_folder = path.join(self.host_folder, type_na)
        self.input_folder = path.join(self.na_folder, 'input')
        self.aa_folder = path.join(self.input_folder, 'allatoms')
        self.heavy_folder = path.join(self.input_folder, 'heavyatoms')
        self.inp_folder = path.join(self.na_folder, 'charmm_inp')
        self.dat_folder = path.join(self.na_folder, 'charmm_dat')
        self.mkcrd_folder = path.join(self.na_folder, 'make_crd')
        self.make_folders()
        self.seq1 = sequences[self.host][self.type_na]['guide']
        self.seq2 = sequences[self.host][self.type_na]['target']
        self.inp_dcd = path.join(self.aa_folder, '{0}.central.dcd'.format(type_na))

    def make_folders(self):
        folders = [self.host_folder, self.na_folder, self.input_folder, self.aa_folder, self.heavy_folder,
                   self.inp_folder, self.dat_folder, self.mkcrd_folder]
        for folder in folders:
            check_dir_exist_and_make(folder)

    def make_crd_input(self, amber=False, firstter=None, lastter=None): 
        if self.type_na == 'arna+arna':
            na = 'arna'
            supplement1 = None
            supplement2 = None 
        elif self.type_na == 'bdna+bdna':
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

        crd2 = path.join(self.mkcrd_folder, '{0}2.crd'.format(na))
        inp2 = Script(path.join(self.mkcrd_folder, '{0}2.inp'.format(na)))
        inp2.write_bomlev()
        inp2.initialize_rtf_prm(amber=amber)
        inp2.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        if supplement2 is not None:
            inp2.write_supplement(supplement2)
        inp2.gen_angle_dihedral()
        inp2.read_pdb(path.join(self.mkcrd_folder, '{0}2.1.pdb'.format(na)))
        inp2.write_crd(crd2)
        inp2.end()

        inp3 = Script(path.join(self.mkcrd_folder, '{0}.inp'.format(self.type_na)))
        inp3.write_bomlev()
        inp3.initialize_rtf_prm(amber=amber)
        inp3.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        if supplement1 is not None:
            inp3.write_supplement(supplement1)
        inp3.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        if supplement2 is not None:
            inp3.write_supplement(supplement2)
        inp3.gen_angle_dihedral()
        inp3.read_crd(crd1, selection='segid strand1', ignore=True)
        inp3.read_crd(crd2, selection='segid strand2', ignore=True)
        inp3.write_crd(path.join(self.mkcrd_folder, '{0}.crd'.format(self.type_na)))
        inp3.end()

    def make_no_h_crd_input(self, amber=False, firstter=None, lastter=None):
        if self.type_na == 'arna+arna':
            supplement1 = None
            supplement2 = None 
        elif self.type_na == 'bdna+bdna':
            supplement1 = get_patch(self.seq1, 1)
            supplement2 = get_patch(self.seq2, 2)

        f_inp = path.join(self.inp_folder, 'write_no_h_crd.inp')
        crd_inp = path.join(self.mkcrd_folder, '{0}.crd'.format(self.type_na))
        crd_out = path.join(self.heavy_folder, '{0}.nohydrogen.crd'.format(self.type_na))

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm(amber=amber)
        inp.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        if supplement1 is not None:
            inp.write_supplement(supplement1)
        inp.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        if supplement2 is not None:
            inp.write_supplement(supplement2)
        inp.gen_angle_dihedral()
        inp.read_crd(crd_inp)
        inp.delete_selection()
        inp.write_crd(crd_out)
        inp.end()

    def make_no_h_dcd_input(self, mass_weighted=True, amber=False, begin=1, frame_num=10000, firstter=None, lastter=None):
        if self.type_na == 'arna+arna':
            supplement1 = None
            supplement2 = None 
        elif self.type_na == 'bdna+bdna':
            supplement1 = get_patch(self.seq1, 1)
            supplement2 = get_patch(self.seq2, 2)

        f_inp = path.join(self.inp_folder, 'write_no_h_dcd.inp')
        crd_inp = path.join(self.mkcrd_folder, '{0}.crd'.format(self.type_na))
        dcd_out = path.join(self.heavy_folder, '{0}.nohydrogen.dcd'.format(self.type_na))

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm(amber=amber)
        inp.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        if supplement1 is not None:
            inp.write_supplement(supplement1)
        inp.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        if supplement2 is not None:
            inp.write_supplement(supplement2)
        inp.gen_angle_dihedral()
        inp.read_crd(crd_inp)
        inp.coor_copy()
        inp.read_traj(self.inp_dcd)
        inp.open_write_traj(dcd_out)
        inp.write_noh_dcd(mass_weighted=mass_weighted, begin=begin, frame_num=frame_num)
        inp.close_unit(21)
        inp.close_unit(30)
        inp.end()

    def make_avg_crd_input(self, selection='all', amber=False, firstter=None, lastter=None):
        if self.type_na == 'arna+arna':
            supplement1 = None
            supplement2 = None 
        elif self.type_na == 'bdna+bdna':
            supplement1 = get_patch(self.seq1, 1)
            supplement2 = get_patch(self.seq2, 2)

        f_inp = path.join(self.inp_folder, 'write_no_h_avg_crd.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.crd'.format(self.type_na))
        dcd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.dcd'.format(self.type_na))
        crd_out = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm(amber=amber)
        inp.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        if supplement1 is not None:
            inp.write_supplement(supplement1)
        inp.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        if supplement2 is not None:
            inp.write_supplement(supplement2)
        inp.gen_angle_dihedral()
        inp.delete_selection()
        inp.read_crd(crd_inp)
        inp.read_traj(dcd_inp)
        inp.calculate_avg(selection=selection)
        inp.coor_copy()
        inp.write_crd(crd_out, comp=True)
        inp.end()


    def fit_dcd_to_avg_input(self, amber=False, begin=1, frame_num=10000, massweighted=False, firstter=None, lastter=None, dcd_out=None):
        if self.type_na == 'arna+arna':
            supplement1 = None
            supplement2 = None 
        elif self.type_na == 'bdna+bdna':
            supplement1 = get_patch(self.seq1, 1)
            supplement2 = get_patch(self.seq2, 2)

        f_inp = path.join(self.inp_folder, 'fit_dcd_to_avg.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        dcd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.dcd'.format(self.type_na))
        if dcd_out is None:
            dcd_out = path.join(self.heavy_folder, '{0}.nohydrogen.fitavg.dcd'.format(self.type_na))

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm(amber=amber)
        inp.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        if supplement1 is not None:
            inp.write_supplement(supplement1)
        inp.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        if supplement2 is not None:
            inp.write_supplement(supplement2)
        inp.gen_angle_dihedral()
        inp.delete_selection()
        inp.read_crd(crd_inp)
        inp.coor_copy()
        inp.read_traj(dcd_inp)
        inp.open_write_traj(dcd_out)
        inp.write_noh_dcd(begin=begin, frame_num=frame_num, mass_weighted=massweighted)
        inp.close_unit(21)
        inp.close_unit(30)
        inp.end()

    def make_convert_avg_to_dcd_input(self, amber=False, firstter=None, lastter=None):
        if self.type_na == 'arna+arna':
            supplement1 = None
            supplement2 = None 
        elif self.type_na == 'bdna+bdna':
            supplement1 = get_patch(self.seq1, 1)
            supplement2 = get_patch(self.seq2, 2)

        f_inp = path.join(self.inp_folder, 'convert_avgcrd_avgdcd.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        dcd_out = path.join(self.heavy_folder, 'aa.avg.dcd')

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm(amber=amber)
        inp.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        if supplement1 is not None:
            inp.write_supplement(supplement1)
        inp.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        if supplement2 is not None:
            inp.write_supplement(supplement2)
        inp.gen_angle_dihedral()
        inp.delete_selection()
        inp.read_crd(crd_inp)
        inp.convert_avg_to_dcd(dcd_out)
        inp.end()

    def make_no_h_crd(self):
        f_inp = path.join(self.inp_folder, 'write_no_h_crd.inp')
        f_dat = path.join(self.dat_folder, 'write_no_h_crd.dat')
        exec_charmm(f_inp, f_dat)

    def make_no_h_dcd(self):
        f_inp = path.join(self.inp_folder, 'write_no_h_dcd.inp')
        f_dat = path.join(self.dat_folder, 'write_no_h_dcd.dat')
        exec_charmm(f_inp, f_dat)

    def make_avg_crd(self):
        f_inp = path.join(self.inp_folder, 'write_no_h_avg_crd.inp')
        f_dat = path.join(self.dat_folder, 'write_no_h_avg_crd.dat')
        exec_charmm(f_inp, f_dat)

    def fit_dcd_to_avg(self):
        f_inp = path.join(self.inp_folder, 'fit_dcd_to_avg.inp')
        f_dat = path.join(self.dat_folder, 'fit_dcd_to_avg.dat')
        exec_charmm(f_inp, f_dat)

    def make_crd(self):
        if self.type_na == 'arna+arna':
            na = 'arna'
        elif self.type_na == 'bdna+bdna':
            na = 'bdna'

        inp1 = path.join(self.mkcrd_folder, '{0}1.inp'.format(na))
        inp1_dat = path.join(self.mkcrd_folder, '{0}1.dat'.format(na))
        exec_charmm(inp1, inp1_dat)

        inp2 = path.join(self.mkcrd_folder, '{0}2.inp'.format(na))
        inp2_dat = path.join(self.mkcrd_folder, '{0}2.dat'.format(na))
        exec_charmm(inp2, inp2_dat)

        inp3 = path.join(self.mkcrd_folder, '{0}.inp'.format(self.type_na))
        inp3_dat = path.join(self.mkcrd_folder, '{0}.dat'.format(self.type_na))
        exec_charmm(inp3, inp3_dat)

    def convert_avg_to_dcd(self):
        f_inp = path.join(self.inp_folder, 'convert_avgcrd_avgdcd.inp')
        f_dat = path.join(self.dat_folder, 'convert_avgcrd_avgdcd.dat')
        exec_charmm(f_inp, f_dat)

    def reset_na2_pdb_resid(self, offset):
        pdb_name = path.join(self.mkcrd_folder, 'bdna2.1.pdb')
        f_backup = path.join(self.mkcrd_folder, 'bdna2.1.backup.pdb')
        copyfile(pdb_name, f_backup)
        print(f'{pdb_name} {f_backup}')
        reader = PDB.PDBReader(pdb_name, skip_header=2, skip_footer=1)
        for atom in reader.atomgroup:
            resid = atom.resid
            atom.set_resid(resid + offset)
        writer = PDB.PDBWriter(pdb_name, reader.atomgroup)
        writer.write_pdb()
        print(f'Reset {pdb_name} resid by offset {offset}!')
        print(f'Check by...\nvim {pdb_name}')


def exec_charmm(f_input, f_output):
    print("charmm< {0} > {1}".format(f_input, f_output))
    check_call(charmm, stdin=open(f_input, 'r'), stdout=open(f_output, 'w+'), shell=True)




