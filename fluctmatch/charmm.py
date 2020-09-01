d_restype = {'nucleic acid': {'A': 'ADE', 'T': 'THY', 'C': 'CYT', 'G': 'GUA', 'U': 'URA'},
             'protein': {}}


class Script:
    def __init__(self, file_name, mode='write', description=None):
        """

        :param file_name:
        :param mode: write, append, modify
        """
        self.file_name = file_name
        if mode == 'write':
            self.f = open(self.file_name, 'w+')
        elif mode == 'append':
            self.f = open(self.file_name, 'a+')
        if description is None:
            self.write_head()
        else:
            self.write_head(author='Yi-Tsao Chen', description=description)

    def write_head(self, author='Yi-Tsao Chen', description=''):
        self.f.write('* Author: {0}\n* {1}\n\n'.format(author, description))

    def write_bomlev(self, bomlev=-2):
        self.f.write('bomlev {0}\n\n'.format(str(bomlev)))

    def initialize_rtf_prm(self, rtfprm_folder='/home/yizaochen/prm_rtf', amber=False):
        self.f.write('set topdir = {0}\n'.format(rtfprm_folder))
        if amber:
            self.f.write('open unit 11 read form name @topdir/parm14sb_all.rtf\n')
            self.f.write('read rtf card unit 11\n')
            self.f.write('close unit 11\n\n')

            self.f.write('open read unit 11 card name @topdir/parm14sb_all.prm\n')
            self.f.write('read param unit 11 flex\n')
            self.f.write('close unit 11\n\n')
        else:
            self.f.write('open unit 11 read form name @topdir/top_all36_prot.rtf\n')
            self.f.write('read rtf card unit 11\n')
            self.f.write('close unit 11\n\n')

            self.f.write('open read unit 11 card name @topdir/par_all36_prot.prm\n')
            self.f.write('read param unit 11 flex\n')
            self.f.write('close unit 11\n\n')

            self.f.write('open unit 11 read form name @topdir/top_all36_na.rtf\n')
            self.f.write('read rtf card unit 11 appe\n')
            self.f.write('close unit 11\n\n')

            self.f.write('open read unit 11 card name @topdir/par_all36_na.prm\n')
            self.f.write('read param unit 11 flex appe\n')
            self.f.write('close unit 11\n\n')

    def write_seq(self, sequence, segid='nucleic', firstter=None, lastter=None, type_molecule='nucleic acid'):
        """

        :param sequence: nucleic or protein sequence, like ATCGGGCT
        :param type_molecule: 'nucleic acid', 'protein'
        :param segid: arna, arna1, arna2, nucleic, protein
        :param firstter: 5ter
        :param lastter: 3ter
        :return:
        """
        self.f.write('read sequ card\n')
        self.f.write('*\n')
        self.f.write(' {0}\n'.format(len(sequence)))
        for res in sequence:
            self.f.write('{0} '.format(d_restype[type_molecule][res]))
        self.f.write('\n\n')
        if (firstter is None) and (lastter is None):
            # self.f.write('generate {0} setup\n\n'.format(segid))
            self.f.write('generate {0} setup first None last None\n\n'.format(segid))
        elif firstter == 'amber_5ter' and lastter == 'amber_3ter':
            self.f.write('generate {0} setup\n\n'.format(segid))
        else:
            self.f.write('generate {0} setup first {1} last {2}\n\n'.format(segid, firstter, lastter))

    def gen_angle_dihedral(self):
        self.f.write('autogenerate angles dihedrals\n\n')

    def delete_selection(self, selection='hydrogen'):
        self.f.write('delete atom sele {0} end\n\n'.format(selection))

    def coor_copy(self):
        self.f.write('coor copy comp\n\n')

    def close_unit(self, unit):
        self.f.write('close unit {0}\n\n'.format(unit))

    def calculate_avg(self, selection='all', begin=1, stop=10000):
        self.f.write('coor dyna firstu 21 begin {0} skip 1 stop {1} sele {2} end\n\n'.format(begin, stop, selection))

    def read_pdb(self, pdbfile):
        self.f.write('open read unit 12 card name {0}\n'.format(pdbfile))
        self.f.write('read coor pdb unit 12\n')
        self.f.write('close unit 12\n\n')

    def read_crd(self, crdfile, selection=None, ignore=True, comp=False):
        self.f.write('open read unit 12 card name {0}\n'.format(crdfile))
        if comp:
            temp = 'read coor comp unit 12'
        else:
            temp = 'read coor unit 12'
        if selection is not None:
            temp += ' sele {0} end'.format(selection)
        if ignore:
            temp += ' ignore'
        temp += '\n'
        self.f.write(temp)
        self.f.write('close unit 12\n\n')

    def write_crd(self, crdfile, comp=False):
        self.f.write('open write unit 12 card name {0}\n'.format(crdfile))
        if comp:
            self.f.write('write coor card comp unit 12\n')
        else:
            self.f.write('write coor card unit 12\n')
        self.f.write('close unit 12\n\n')

    def read_traj(self, trajfile):
        self.f.write('open unform read unit 21 name {0}\n'.format(trajfile))
        self.f.write('traj firstu 21 nunit 1 skip 1\n\n')

    def open_write_traj(self, trajfile):
        self.f.write('open unform write unit 30 name {0}\n\n'.format(trajfile))

    def write_noh_dcd(self, mass_weighted=True, begin=1, frame_num=10000):
        self.f.write('merge coor firstu 21 nunit 1 outputu 30 -\n')
        self.f.write('begin {0} skip 1 nfile {1} -\n'.format(begin, frame_num))
        self.f.write('sele .not. hydrogen end -\n')
        if mass_weighted:
            self.f.write('orie mass sele .not. hydrogen end\n\n')
        else:
            self.f.write('orie sele .not. hydrogen end\n\n')

    def write_quasi(self, vect=True):
        self.f.write('set nmod ?natom\n')
        self.f.write('calc nmod = @nmod * 3\n\n')
        self.f.write('vibran nmod @nmod\n')
        self.f.write('quasi temp 310 firstu 21 nunit 1 begin 1 skip 1\n')
        if vect:
            self.f.write('print norm vect\n\n')
        else:
            self.f.write('\n')

    def write_covar_mat_crd(self, crdfile):
        self.f.write('open write unit 15 card name {0}\n'.format(crdfile))
        self.f.write('coor cova firstu 21 nunit 1 matrix unit 15\n\n')

    def write_quasi_mode_crd(self, crdfile):
        self.f.write('open write unit 15 card name {0}\n'.format(crdfile))
        self.f.write('write norm card mode 1 thru @nmod unit 15\n')
        self.f.write('close unit 15\n\n')

    def write_quasi_single_mode_crd(self, modeid, crdfile):
        self.f.write('set modeid {0}\n'.format(modeid))
        self.f.write('calc modeid1 = @nmod - @modeid + 1\n')
        self.f.write('open write unit 15 card name {0}\n'.format(crdfile))
        self.f.write('write norm card mode @modeid1 unit 15\n')
        self.f.write('close unit 15\n\n')

    def write_quasi_single_mode_crd_fm(self, modeid, crdfile, lastsixmode=False):
        self.f.write('set modeid {0}\n'.format(modeid))
        if lastsixmode:
            self.f.write('calc modeid1 = @modeid\n')
        else:
            self.f.write('calc modeid1 = @modeid + 6\n')
        self.f.write('open write unit 15 card name {0}\n'.format(crdfile))
        self.f.write('write norm card mode @modeid1 unit 15\n')
        self.f.write('close unit 15\n\n')

    def read_normal_mode(self, crdfile):
        self.f.write('set nmod ?natom\n')
        self.f.write('calc nmod = @nmod * 3\n')
        self.f.write('vibran nmod @nmod\n\n')
        self.f.write('open read unit 15 card name {0}\n'.format(crdfile))
        self.f.write('read norm card unit 15 mode 1 thru @nmod\n')
        self.f.write('close unit 15\n\n')

    def write_quasi_mode_dcd(self, dcdfile, modeid=1, phase='3.6', ncyc=10, temp=310, tfre=1):
        self.f.write('set modeid {0}\n'.format(modeid))
        self.f.write('calc modeid1 = @nmod - @modeid + 1\n')
        self.f.write('open unit 15 write unform name {0}\n'.format(dcdfile))
        self.f.write('write traj mode @modeid1 step 0.1 phas {0} ncyc {1} temp {2} tfre {3} unit 15\n'.format(phase,
                     ncyc, temp, tfre))
        self.f.write('close unit 15\n\n')

    def write_quasi_mode_supe_dcd(self, dcdfile, modeid=1):
        self.f.write('set modeid {0}\n'.format(modeid))
        self.f.write('calc modeid1 = @nmod - @modeid + 1\n')
        self.f.write('open unit 15 write unform name {0}\n'.format(dcdfile))
        if modeid == 1:
            self.f.write('write traj mode @modeid1 step 0.1 phas 3.6 ncyc 10 temp 310 unit 15\n')
        else:
            self.f.write('write traj mode @modeid1 thru @nmod supe step 0.1 phas 3.6 ncyc 10 temp 310 unit 15\n')
        self.f.write('close unit 15\n\n')

    def write_proj_mode_crd(self, crdfile, modeid=1):
        self.f.write('set modeid {0}\n'.format(modeid))
        self.f.write('calc modeid1 = @nmod - @modeid + 1\n')
        self.f.write('open write unit 16 card name {0}\n'.format(crdfile))
        self.f.write('proj traj nunit 1 firstu 21 mode @modeid1 write wunit 16 card\n')
        self.f.write('close unit 16\n\n')

    def make_diff_to_avg(self, dcdin, dcdout, n_frame):
        self.f.write('set framenum {0}\n'.format(n_frame))
        self.f.write('calc criteria = @framenum + 1\n\n')
        self.f.write('open unform read unit 21 name {0}\n'.format(dcdin))
        self.f.write('open unform write unit 30 name {0}\n'.format(dcdout))
        self.f.write('traj firstu 21 iwrite 30 nfile @framenum\n\n')
        self.f.write('set i 1\n')
        self.f.write('label loop\n')
        self.f.write('traj read\n')
        self.f.write('coor diff\n')
        self.f.write('traj write\n')
        self.f.write('incr i by 1\n')
        self.f.write('if @i lt @criteria goto loop\n\n')

    def write_dcd_by_selected_frame(self, frame_id, dcd_in, dcd_out):
        self.f.write('set selectedframe {0}\n'.format(frame_id))
        self.f.write('open unform read unit 21 name {0}\n'.format(dcd_in))
        self.f.write('traj firstu 21 nunit 1 begin @selectedframe skip 1\n\n')
        self.f.write('set framenum 1\n')
        self.f.write('traj read\n')
        self.f.write('open unform write unit 30 name {0}\n'.format(dcd_out))
        self.f.write('traj iwrite 30 nfile @framenum\n')
        self.f.write('traj write\n\n')

    def proj_all_modes_to_crd(self, crd_out):
        self.f.write('open write unit 16 card name {0}\n'.format(crd_out))
        self.f.write('set modeid 1\n')
        self.f.write('label loop\n')
        self.f.write('calc modeid1 = @nmod - @modeid + 1\n')
        self.f.write('proj traj nunit 1 firstu 21 mode @modeid1 write wunit 16 card\n')
        self.f.write('incr modeid BY 1\n')
        self.f.write('if @modeid lt 1267 goto loop\n')
        self.f.write('close unit 16\n\n')

    def print_eigenvector(self, modeid=1):
        self.f.write('set modeid {0}\n'.format(modeid))
        self.f.write('calc modeid1 = @nmod - @modeid + 1\n')
        self.f.write('print norm mode @modeid1 vect\n\n')

    def set_all_mass_same(self):
        self.f.write('scalar mass set 1 sele all end\n\n')

    def quasi_end(self):
        self.f.write('end\n\n')

    def write_supplement(self, supplement):
        self.f.write(supplement)
        self.f.write('\n')

    def convert_avg_to_dcd(self, dcd_out):
        self.f.write('coor copy comp\n')
        self.f.write('set framenum 1\n')
        self.f.write('open unform write unit 30 name {0}\n'.format(dcd_out))
        self.f.write('traj iwrite 30 nfile @framenum\n')
        self.f.write('traj write')

    def end(self):
        self.f.write('stop')
        self.f.close()

