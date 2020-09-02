from os import path, mkdir
from subprocess import check_call, check_output
from charmm import Script
from numpy import genfromtxt
import pandas as pd
import matplotlib.pyplot as plt
import re
import util

root = '/Users/yizao/IdeaProjects/python/connect_macro_micro'
charmm = '/Users/yizao/charmm/exec/osx/charmm'
d_sequence = {'pnas': {'arna+arna': {1: 'CAAUGGAGUA', 2: 'UACUCCAUUG'},
                       'bdna+bdna': {1: 'CAATGGAGTA', 2: 'TACTCCATTG'}},
              'pnas_amber': {'arna+arna': {1: 'CAAUGGAGUA', 2: 'UACUCCAUUG'},
                             'bdna+bdna': {1: 'CAATGGAGTA', 2: 'TACTCCATTG'}},
              'pnas_clean': {'arna+arna': {1: 'CAAUGGAGUA', 2: 'UACUCCAUUG'},
                             'bdna+bdna': {1: 'CAATGGAGTA', 2: 'TACTCCATTG'}},
              'pnas_amber_clean': {'arna+arna': {1: 'CAAUGGAGUA', 2: 'UACUCCAUUG'},
                                   'bdna+bdna': {1: 'CAATGGAGTA', 2: 'TACTCCATTG'}},
              'pnas_amber_20ns': {'arna+arna': {1: 'CAAUGGAGUA', 2: 'UACUCCAUUG'},
                                   'bdna+bdna': {1: 'CAATGGAGTA', 2: 'TACTCCATTG'}},
              'pnas_amber_16mer': {'arna+arna': {1: 'GCGCAAUGGAGUACGC', 2: 'GCGUACUCCAUUGCGC'},
                                   'bdna+bdna': {1: 'GCGCAATGGAGTACGC', 2: 'GCGTACTCCATTGCGC'}}
              }


class QuasiAgent:
    def __init__(self, host, type_na):
        self.host = host
        self.type_na = type_na
        self.host_folder = path.join(root, host)
        self.na_folder = path.join(self.host_folder, type_na)
        self.input_folder = path.join(self.na_folder, 'input')
        self.aa_folder = path.join(self.input_folder, 'allatoms')
        self.heavy_folder = path.join(self.input_folder, 'heavyatoms')
        self.inp_folder = path.join(self.na_folder, 'charmm_inp')
        self.dat_folder = path.join(self.na_folder, 'charmm_dat')
        self.mkcrd_folder = path.join(self.na_folder, 'make_crd')
        self.nm_folder = path.join(self.na_folder, 'nm')
        self.nm_single_folder = path.join(self.na_folder, 'nm_single_mode')
        self.mode_traj_folder = path.join(self.na_folder, 'mode_traj')
        self.covar_folder = path.join(self.na_folder, 'covar_mat')
        self.freq_folder = path.join(self.na_folder, 'frequency')
        self.coef_folder = path.join(self.na_folder, 'proj_coef')
        self.diff_folder = path.join(self.na_folder, 'diff_crd')
        self.eigenvector_folder = path.join(self.na_folder, 'eigenvectors')
        self.make_folders()
        self.inp_dcd = path.join(self.aa_folder, '{0}.central.dcd'.format(type_na))
        self.check_dcd()
        self.template_folder = path.join(root, 'charmm_template')
        if re.match('pnas_amber_clean', self.host):
            seq_host = 'pnas_amber_clean'
            self.seq1 = d_sequence[seq_host][self.type_na][1]
            self.seq2 = d_sequence[seq_host][self.type_na][2]
        elif re.match('pnas_amber_16mer', self.host):
            seq_host = 'pnas_amber_16mer'
            self.seq1 = d_sequence[seq_host][self.type_na][1]
            self.seq2 = d_sequence[seq_host][self.type_na][2]
        else:
            self.seq1 = d_sequence[self.host][self.type_na][1]
            self.seq2 = d_sequence[self.host][self.type_na][2]

    def make_folders(self):
        folders = [self.host_folder, self.na_folder, self.input_folder, self.aa_folder, self.heavy_folder,
                   self.inp_folder, self.dat_folder, self.mkcrd_folder, self.nm_folder, self.mode_traj_folder,
                   self.covar_folder, self.freq_folder, self.coef_folder, self.diff_folder, self.eigenvector_folder,
                   self.nm_single_folder]
        for folder in folders:
            check_dir_exist_and_make(folder)

    def write_freq_file(self, mode):
        f_in = path.join(self.dat_folder, 'quasi_mode.dat')
        f_out = path.join(self.freq_folder, 'mode.{0}.txt'.format(mode))
        cmd = 'grep -A 256 FREQUENCIES {0} > {1}'.format(f_in, f_out)
        check_call(cmd, shell=True)

    def check_dcd(self):
        if not path.exists(self.inp_dcd):
            print("Convert XTC to DCD....")
            xtc = path.join(self.aa_folder, '{0}.central.xtc'.format(self.type_na))
            util.xtc2dcd(xtc, self.inp_dcd)
        else:
            print("{0} already exists".format(self.inp_dcd))

    def make_crd_input(self, amber=False, supplement1=None, supplement2=None, firstter=None, lastter=None, offset=-3): 
        if self.type_na == 'arna+arna':
            na = 'arna'
        elif self.type_na == 'bdna+bdna':
            na = 'bdna'

        pdb1_1 = path.join(self.mkcrd_folder, '{0}1.1.pdb'.format(na))
        pdb2_1 = path.join(self.mkcrd_folder, '{0}2.1.pdb'.format(na))
        pdb1_2 = path.join(self.mkcrd_folder, '{0}1.2.pdb'.format(na))
        pdb2_2 = path.join(self.mkcrd_folder, '{0}2.2.pdb'.format(na))
        util.change_resid(pdb1_1, offset, pdb1_2)
        util.change_resid(pdb2_1, offset, pdb2_2)

        crd1 = path.join(self.mkcrd_folder, '{0}1.crd'.format(na))
        inp1 = Script(path.join(self.mkcrd_folder, '{0}1.inp'.format(na)))
        inp1.write_bomlev()
        inp1.initialize_rtf_prm(amber=amber)
        inp1.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        if supplement1 is not None:
            inp1.write_supplement(supplement1)
        inp1.gen_angle_dihedral()
        inp1.read_pdb(path.join(self.mkcrd_folder, '{0}1.2.pdb'.format(na)))
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
        inp2.read_pdb(path.join(self.mkcrd_folder, '{0}2.2.pdb'.format(na)))
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

    def make_no_h_crd_input(self, amber=False, supplement1=None, supplement2=None, firstter=None, lastter=None):
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

    def make_no_h_dcd_input(self, mass_weighted=True, amber=False, supplement1=None, supplement2=None,
                            begin=1, frame_num=10000, firstter=None, lastter=None):
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

    def make_avg_crd_input(self, selection='all', amber=False, supplement1=None, supplement2=None, firstter=None, lastter=None):
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

    def make_get_eigenvector_input(self, modeid=1, firstter=None, lastter=None):
        f_inp = path.join(self.inp_folder, 'get_eigenvectors.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        mode_inp = path.join(self.nm_folder, 'quasi.nmvector.crd')

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm()
        inp.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        inp.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        inp.gen_angle_dihedral()
        inp.delete_selection()
        inp.read_crd(crd_inp)
        inp.read_normal_mode(mode_inp)
        inp.print_eigenvector(modeid)
        inp.quasi_end()
        inp.end()

    def make_get_eigenvector_input_v1(self, modeid=1, amber=False, supplement1=None, supplement2=None, firstter=None, lastter=None):
        f_inp = path.join(self.inp_folder, 'get_eigenvectors.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        mode_inp = path.join(self.nm_folder, 'quasi.nmvector.crd')
        crd_out = path.join(self.nm_single_folder, 'mode.{0}.crd'.format(modeid))

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
        inp.read_normal_mode(mode_inp)
        inp.write_quasi_single_mode_crd(modeid, crd_out)
        inp.quasi_end()
        inp.end()

    def fit_dcd_to_avg_input(self, amber=False, supplement1=None, supplement2=None, begin=1, frame_num=10000,
                             massweighted=False, firstter=None, lastter=None, dcd_out=None):
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

    def make_quasi_input(self, mass_all_same=False, amber=False, supplement1=None, supplement2=None,
                         begin=1, stop=10000, firstter=None, lastter=None):
        f_inp = path.join(self.inp_folder, 'quasi.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        dcd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.fitavg.dcd'.format(self.type_na))
        crd_out = path.join(self.nm_folder, 'quasi.nmvector.crd')

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
        if mass_all_same:
            inp.set_all_mass_same()
        inp.read_crd(crd_inp)
        inp.read_traj(dcd_inp)
        inp.calculate_avg(begin=begin, stop=stop)
        inp.coor_copy()
        inp.write_quasi()
        inp.write_quasi_mode_crd(crd_out)
        inp.quasi_end()
        inp.close_unit(21)
        inp.end()

    def make_quasi_mode_input(self, mode=1, firstter=None, lastter=None):
        f_inp = path.join(self.inp_folder, 'quasi_mode.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        dcd_inp = path.join(self.mode_traj_folder, 'mode.{0}.dcd'.format(str(mode)))

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm()
        inp.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        inp.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        inp.gen_angle_dihedral()
        inp.delete_selection()
        inp.read_crd(crd_inp)
        inp.read_traj(dcd_inp)
        inp.calculate_avg()
        inp.coor_copy()
        inp.write_quasi(vect=False)
        inp.quasi_end()
        inp.close_unit(21)
        inp.end()

    def make_gen_diff_crd(self, mode=1, n_frame=1000, firstter=None, lastter=None):
        f_inp = path.join(self.inp_folder, 'gen_diff_crd.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        dcd_inp = path.join(self.mode_traj_folder, 'mode.{0}.dcd'.format(str(mode)))
        dcd_out = path.join(self.diff_folder, 'mode.{0}.dcd'.format(str(mode)))

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm()
        inp.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        inp.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        inp.gen_angle_dihedral()
        inp.delete_selection()
        inp.read_crd(crd_inp)
        inp.coor_copy()
        inp.make_diff_to_avg(dcd_inp, dcd_out, n_frame)
        inp.end()

    def make_proj_mode_input(self, mode=1):
        f_inp = path.join(self.inp_folder, 'proj_mode.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        mode_inp = path.join(self.nm_folder, 'quasi.nmvector.crd')
        dcd_inp = path.join(self.diff_folder, 'mode.{0}.dcd'.format(str(mode)))
        crd_out = path.join(self.coef_folder, 'proj.mode.{0}.crd'.format(str(mode)))

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm()
        inp.write_seq(self.seq1, segid='strand1')
        inp.write_seq(self.seq2, segid='strand2')
        inp.gen_angle_dihedral()
        inp.delete_selection()
        inp.read_crd(crd_inp)
        # inp.coor_copy()
        inp.read_traj(dcd_inp)
        inp.read_normal_mode(mode_inp)
        inp.write_proj_mode_crd(crd_out, mode)
        inp.quasi_end()
        inp.close_unit(21)
        inp.end()

    def make_proj_md_to_mode_input(self, mode):
        f_inp = path.join(self.inp_folder, 'proj_md_to_mode.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        mode_inp = path.join(self.nm_folder, 'quasi.nmvector.crd')
        dcd_inp = path.join(self.diff_folder, 'mode.0.dcd')
        crd_out = path.join(self.coef_folder, 'proj.md.onto.mode.{0}.crd'.format(str(mode)))

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm()
        inp.write_seq(self.seq1, segid='strand1')
        inp.write_seq(self.seq2, segid='strand2')
        inp.gen_angle_dihedral()
        inp.delete_selection()
        inp.read_crd(crd_inp)
        # inp.coor_copy()
        inp.read_traj(dcd_inp)
        inp.read_normal_mode(mode_inp)
        inp.write_proj_mode_crd(crd_out, mode)
        inp.quasi_end()
        inp.close_unit(21)
        inp.end()

    def make_covariance_mat_input(self, modeid=0):
        f_inp = path.join(self.inp_folder, 'write_covar_mat.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        dcd_inp = path.join(self.mode_traj_folder, 'mode.{0}.dcd'.format(str(modeid)))
        crd_out = path.join(self.covar_folder, 'covar.{0}.crd'.format(str(modeid)))

        inp = Script(f_inp)
        inp.write_bomlev()
        inp.initialize_rtf_prm()
        inp.write_seq(self.seq1, segid='strand1')
        inp.write_seq(self.seq2, segid='strand2')
        inp.gen_angle_dihedral()
        inp.delete_selection()
        inp.read_crd(crd_inp)
        inp.read_traj(dcd_inp)
        inp.write_covar_mat_crd(crd_out)
        inp.close_unit(15)
        inp.close_unit(21)
        inp.end()

    def make_gen_mode_input(self, modeid, phase='3.6', ncyc=10, temp=310, tfre=1, amber=False,
                            supplement1=None, supplement2=None, firstter=None, lastter=None):
        f_inp = path.join(self.inp_folder, 'gen_mode.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        mode_inp = path.join(self.nm_folder, 'quasi.nmvector.crd')
        dcd_out = path.join(self.mode_traj_folder, 'mode.{0}.dcd'.format(str(modeid)))

        inp = Script(f_inp, mode='write', description='generate mode {0}'.format(modeid))
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
        inp.read_normal_mode(mode_inp)
        inp.write_quasi_mode_dcd(dcd_out, modeid, phase=phase, ncyc=ncyc, temp=temp, tfre=tfre)
        inp.quasi_end()
        inp.end()

    def make_convert_avg_to_dcd_input(self, amber=False, supplement1=None, supplement2=None, firstter=None, lastter=None):
        f_inp = path.join(self.inp_folder, 'convert_avgcrd_avgdcd.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        dcd_out = path.join(self.mode_traj_folder, 'aa.avg.dcd')

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

    def make_gen_mode_superimpose_input(self, modeid):
        f_inp = path.join(self.inp_folder, 'gen_mode.inp')
        crd_inp = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        mode_inp = path.join(self.nm_folder, 'quasi.nmvector.crd')
        dcd_out = path.join(self.mode_traj_folder, 'mode.1to{0}.dcd'.format(str(modeid)))

        inp = Script(f_inp, mode='write', description='generate mode {0}'.format(modeid))
        inp.write_bomlev()
        inp.initialize_rtf_prm()
        inp.write_seq(self.seq1, segid='strand1')
        inp.write_seq(self.seq2, segid='strand2')
        inp.gen_angle_dihedral()
        inp.delete_selection()
        inp.read_crd(crd_inp)
        inp.read_normal_mode(mode_inp)
        inp.write_quasi_mode_supe_dcd(dcd_out, modeid)
        inp.quasi_end()
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

    def do_quasi(self):
        f_inp = path.join(self.inp_folder, 'quasi.inp')
        f_dat = path.join(self.dat_folder, 'quasi.dat')
        exec_charmm(f_inp, f_dat)

    def do_quasi_mode(self):
        f_inp = path.join(self.inp_folder, 'quasi_mode.inp')
        f_dat = path.join(self.dat_folder, 'quasi_mode.dat')
        exec_charmm(f_inp, f_dat)

    def gen_mode(self):
        f_inp = path.join(self.inp_folder, 'gen_mode.inp')
        f_dat = path.join(self.dat_folder, 'gen_mode.dat')
        exec_charmm(f_inp, f_dat)

    def gen_diff_crd(self):
        f_inp = path.join(self.inp_folder, 'gen_diff_crd.inp')
        f_dat = path.join(self.dat_folder, 'gen_diff_crd.dat')
        exec_charmm(f_inp, f_dat)

    def make_covar_mat(self):
        f_inp = path.join(self.inp_folder, 'write_covar_mat.inp')
        f_dat = path.join(self.dat_folder, 'write_covar_mat.dat')
        exec_charmm(f_inp, f_dat)

    def proj_mode(self):
        f_inp = path.join(self.inp_folder, 'proj_mode.inp')
        f_dat = path.join(self.dat_folder, 'proj_mode.dat')
        exec_charmm(f_inp, f_dat)

    def proj_md_to_mode(self):
        f_inp = path.join(self.inp_folder, 'proj_md_to_mode.inp')
        f_dat = path.join(self.dat_folder, 'proj_md_to_mode.dat')
        exec_charmm(f_inp, f_dat)

    def get_eigenvector(self):
        f_inp = path.join(self.inp_folder, 'get_eigenvectors.inp')
        f_dat = path.join(self.dat_folder, 'get_eigenvectors.dat')
        exec_charmm(f_inp, f_dat)

    def convert_avg_to_dcd(self):
        f_inp = path.join(self.inp_folder, 'convert_avgcrd_avgdcd.inp')
        f_dat = path.join(self.dat_folder, 'convert_avgcrd_avgdcd.dat')
        exec_charmm(f_inp, f_dat)

    def show_mode(self, mode, superimpose=False):
        crd = path.join(self.heavy_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        if superimpose:
            dcd = path.join(self.mode_traj_folder, 'mode.1to{0}.dcd'.format(str(mode)))
        else:
            dcd = path.join(self.mode_traj_folder, 'mode.{0}.dcd'.format(str(mode)))
        return crd, dcd

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

    def make_eigenvector(self, n_atoms, modeid):
        f_in = path.join(self.dat_folder, 'get_eigenvectors.dat')
        f_out = path.join(self.eigenvector_folder, 'mode.{0}.crd'.format(modeid))
        cmd = 'grep -A {0} \'EIGENVECTOR\' {1} | sed \'s/^            //\' > {2}'.format(n_atoms, f_in, f_out)
        print(cmd)
        check_call(cmd, shell=True)

    def make_eigenvector_v1(self, modeid):
        f_in = path.join(self.nm_single_folder, 'mode.{0}.crd'.format(modeid))
        f_out = path.join(self.nm_single_folder, 'mode.{0}.txt'.format(modeid))
        f = open(f_in, 'r')
        lines = f.readlines()
        f.close()
        f = open(f_out, 'w')
        for line in lines[76:]:
            f.write('{0}'.format(line))
        f.close()


def exec_charmm(f_input, f_output):
    print("charmm< {0} > {1}".format(f_input, f_output))
    check_call(charmm, stdin=open(f_input, 'r'), stdout=open(f_output, 'w+'), shell=True)


def check_dir_exist_and_make(file_path):
    if path.exists(file_path):
        print("{0} exists".format(file_path))
    else:
        print("mkdir {0}".format(file_path))
        mkdir(file_path)


def make_covar_df(host, type_na, modelist, f_output=None):
    from numpy.linalg import norm
    import pandas as pd
    data_folder = path.join(root, host, type_na, 'covar_mat')
    d = {'mode': list(), 'norm': list()}
    for mode in modelist:
        f = path.join(data_folder, 'covar.{0}.crd'.format(str(mode)))
        mat = genfromtxt(f)
        norm_mat = norm(mat)
        d['mode'].append(mode)
        d['norm'].append(norm_mat)
    df = pd.DataFrame(d)
    if f_output is not None:
        df.to_csv(f_output, index=False)
    return df


def get_covar_mat(host, type_na, mode=0):
    data_folder = path.join(root, host, type_na, 'covar_mat')
    f = path.join(data_folder, 'covar.{0}.crd'.format(str(mode)))
    return genfromtxt(f)


def get_last_frequency(f_input):
    cmd = 'grep -A 256 FREQUENCIES {0} | grep \' 1272\''.format(f_input)
    grep_out = check_output(cmd, shell=True)
    return float(grep_out.split()[-1])


def get_raw_frequency_dict(host, type_na):
    f_in = path.join(root, host, type_na, 'frequency', 'mode.0.txt')
    f = open(f_in, 'r')
    lines = f.readlines()
    f.close()
    lines = lines[2:]
    d = dict()
    for line in lines:
        temp = line.split()
        for key, value in zip(temp[::2],temp[1::2]):
            mode = 1273 - int(key)
            d[mode] = float(value)
    return d


def get_frequency_df(host, type_na, modelist):
    d = {'mode': list(), 'raw': list(), 'frequency': list()}
    d_raw = get_raw_frequency_dict(host, type_na)
    for mode in modelist:
        f_in = path.join(root, host, type_na, 'frequency', 'mode.{0}.txt'.format(mode))
        freq = get_last_frequency(f_in)
        d['mode'].append(mode)
        d['raw'].append(d_raw[mode])
        d['frequency'].append(freq)
    df = pd.DataFrame(d)
    f_output = path.join(root, host, type_na, 'frequency', 'single_mode_freq.csv')
    df.to_csv(f_output, index=False)
    return df[['mode', 'raw', 'frequency']]


def plot_raw_frequency(host, type_na):
    f_in = path.join(root, host, type_na, 'frequency', 'single_mode_freq.csv')
    df = pd.read_csv(f_in)
    x = df['mode'].tolist()
    raw = df['raw'].tolist()
    freq = df['frequency'].tolist()
    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(10, 4))

    ax1.plot(x, raw, label='Raw')
    ax1.plot(x, freq, 'ro', label='Quasi-Mode', markersize=0.5)
    ax1.set_xlabel('mode', fontsize=14)
    ax1.set_ylabel('frequency', fontsize=14)
    ax1.legend(loc='best')

    ax2.plot(x, raw, label='Raw')
    ax2.plot(x, freq, 'ro', label='Quasi-Mode')
    ax2.set_xlabel('mode', fontsize=14)
    ax2.set_ylabel('frequency', fontsize=14)
    ax2.set_xlim(0, 10)
    ax2.set_ylim(1, 10)
    ax2.legend(loc='best')
    plt.show()


