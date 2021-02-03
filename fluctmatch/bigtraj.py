from os import path, system, remove
from shutil import copyfile, make_archive
from fluctmatch.avg_dcd_noh import AvgcrddcdAgent, exec_charmm
from fluctmatch.miscell import check_dir_exist_and_make, get_patch
from fluctmatch.sequence import sequences
from fluctmatch.charmm import Script
from fluctmatch.enm import ENMAgent
from fluctmatch.rtf import RTF
from fluctmatch.ic_str import ICSTR

class BigTrajAgent(AvgcrddcdAgent):

    start_time = 0
    end_time = 5000 # 5000 ns
    interval_time = 1000 # unit: ns
    gmx = '/usr/bin/gmx'

    multiscale_bigtraj_folder = '/home/yizaochen/bigtraj_fluctmatch'

    def __init__(self, host, type_na, allsys_folder, bigtraj_folder, simu_folder):
        super().__init__(host, type_na, allsys_folder)
        self.allsys_folder = allsys_folder
        self.bigtraj_folder = bigtraj_folder
        self.simu_folder = simu_folder

        self.time_list, self.mdnum_list = self.get_time_list()
        self.d_smallagents = self.get_all_small_agents()

        self.d_pairs = None

    def get_time_list(self):
        middle_interval = int(self.interval_time/2)
        time_list = list()
        mdnum_list = list()
        mdnum1 = 1
        for time1 in range(self.start_time, self.end_time, middle_interval):
            time2 = time1 + self.interval_time
            if time2 <= self.end_time:
                time_list.append((time1, time2))
                mdnum_list.append((mdnum1, mdnum1+9))
            mdnum1 += 5
        return time_list, mdnum_list

    def get_all_small_agents(self):
        d_smallagents = dict()
        for time1, time2 in self.time_list:
            time_label = f'{time1}_{time2}'
            d_smallagents[(time1,time2)] = SmallTrajAgent(self.bigtraj_folder, self.host, self.type_na, time_label)
        return d_smallagents

    def initialize_all_small_folders(self):
        for time1, time2 in self.time_list:
            self.d_smallagents[(time1,time2)].initialize_folders()

    def copy_refcrd_to_smallfolders(self):
        refcrd = path.join(self.heavy_folder, '{0}.nohydrogen.crd'.format(self.type_na))
        for time1, time2 in self.time_list:
            self.d_smallagents[(time1,time2)].get_refcrd(refcrd)

    def concatenate_xtc_by_gmx(self):
        for timezip, mdnum_zip in zip(self.time_list, self.mdnum_list):
            time1, time2 = timezip
            mdnum1, mdnum2 = mdnum_zip
            self.d_smallagents[(time1,time2)].concatenate_trajectory(self.gmx, self.simu_folder, self.type_na, mdnum1, mdnum2)
    
    def remove_all_redudant_xtc_dcd(self):
        for time1, time2 in self.time_list:
            self.d_smallagents[(time1,time2)].remove_redundant_trajectories()

    def convert_xtc_to_dcd_by_vmd(self):
        for time1, time2 in self.time_list:
            refpdb = path.join(self.allsys_folder, self.host, self.type_na, 'input', 'allatoms', f'{self.type_na}.npt4.all.pdb')
            self.d_smallagents[(time1,time2)].vmd_xtc_to_dcd(refpdb)

    def check_vmd_dcd_status(self):
        for time1, time2 in self.time_list:
            self.d_smallagents[(time1,time2)].check_size('dcd_in')

    def remove_hydrogen_by_charmm(self):
        begin_frame = 1
        frame_num = 10000
        crd_inp = path.join(self.allsys_folder, self.host, self.type_na, 'make_crd', f'{self.type_na}.crd')
        for time1, time2 in self.time_list:
            self.d_smallagents[(time1,time2)].make_no_h_dcd_input(crd_inp, amber=True, begin=begin_frame, frame_num=frame_num, firstter='amber_5ter', lastter='amber_3ter')
            self.d_smallagents[(time1,time2)].make_no_h_dcd()

    def check_nohydrogen_dcd_status(self):
        for time1, time2 in self.time_list:
            self.d_smallagents[(time1,time2)].check_size('dcd_out')

    def set_required_dictionaries(self):
        for time1, time2 in self.time_list:
            self.d_smallagents[(time1,time2)].set_mda_universe()
            self.d_smallagents[(time1,time2)].set_required_d()

    def get_all_pairs(self, cutoff):
        self.d_pairs = dict()
        for time1, time2 in self.time_list:
            pairs = self.d_smallagents[(time1,time2)].get_all_pairs(cutoff)
            n_atoms = len(self.d_smallagents[(time1,time2)].atomid_map)
            n_pairs = len(pairs)
            tl = self.d_smallagents[(time1,time2)].time_label
            print(f'{tl}: {n_atoms} beads and {n_pairs} pairs')
            self.d_pairs[(time1,time2)] = pairs

    def make_all_rtfs(self, cutoff):
        for time1, time2 in self.time_list:
            agent = self.d_smallagents[(time1,time2)]
            pairs = self.d_pairs[(time1,time2)]
            rtf_agent = RTF(self.host, self.type_na, agent.crd, pairs, agent.mass_map)
            f_out = path.join(agent.rtficstr_folder, f'na_enm_{cutoff:.2f}.rtf')
            rtf_agent.write_rtf(f_out)
            print(f'Generate {f_out}')

    def make_all_icstrs(self, cutoff):
        for time1, time2 in self.time_list:
            agent = self.d_smallagents[(time1,time2)]
            pairs = self.d_pairs[(time1,time2)]
            f_out = path.join(agent.rtficstr_folder, f'na_enm_{cutoff:.2f}.str')
            icstr_agent = ICSTR(self.host, self.type_na, agent.crd, pairs, agent.atomid_map)
            icstr_agent.write_ic_str(f_out)
            print(f'Generate {f_out}')

    def make_all_enm_crds(self):
        for time1, time2 in self.time_list:
            agent = self.d_smallagents[(time1,time2)]
            agent.write_make_enm_crd_input(firstter='amber_5ter', lastter='amber_3ter')
            agent.make_enm_crd()

    def archive_host_folder(self):
        output_name = path.join(self.bigtraj_folder, f'{self.host}')
        target_name = path.join(self.bigtraj_folder, self.host)
        make_archive(output_name, 'zip', target_name)
        print(f'Archive {target_name} into {output_name}.zip')

    def update_to_multiscale(self):
        old_f = path.join(self.bigtraj_folder, f'{self.host}.zip')
        new_f = path.join(self.multiscale_bigtraj_folder, f'{self.host}.zip')
        cmd = f'scp {old_f} yizaochen@multiscale:{new_f}'
        print(cmd)

class SmallTrajAgent(ENMAgent):
    scratchroot = '/scratch'

    def __init__(self, root, host, type_na, time_label):
        self.host = host
        self.type_na = type_na
        self.time_label = time_label

        self.host_folder = path.join(root, host)
        self.na_folder = path.join(self.host_folder, type_na)
        self.time_folder = path.join(self.na_folder, time_label) 
        self.input_folder = path.join(self.time_folder, 'input')
        self.charmminp_folder = path.join(self.time_folder, 'charmm_inp')
        self.charmmdat_folder = path.join(self.time_folder, 'charmm_dat')
        #self.mode_traj_folder = path.join(self.time_folder, 'mode_traj')
        self.ic_folder = path.join(self.time_folder, 'ic')
        self.mat_folder = path.join(self.time_folder, 'ic_fluct_mat')
        self.rtficstr_folder = path.join(self.time_folder, 'rtf_ic_str')
        self.datafolder = path.join(self.time_folder, 'data')
        self.backupfolder = path.join(self.datafolder, 'backup')
        #self.initialize_folders()

        # Scratch folder
        self.scr_host = path.join(self.scratchroot, host)
        self.scr_na = path.join(self.scr_host, type_na)
        self.scratchfolder = path.join(self.scr_na, self.time_label)

        self.seq1 = sequences[self.host][self.type_na]['guide']
        self.seq2 = sequences[self.host][self.type_na]['target']

        self.enmprm = path.join(self.datafolder, 'na_enm.prm')
        #self.mode0dcd = path.join(self.mode_traj_folder, 'mode.0.dcd')
        #self.crd = path.join(self.input_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
        self.crd = path.join(self.input_folder, '{0}.nohydrogen.crd'.format(self.type_na))
        self.xtc_in = path.join(self.input_folder, f'{time_label}.xtc')
        self.dcd_in = path.join(self.input_folder, f'{time_label}.dcd')
        self.dcd_out = path.join(self.input_folder, f'{time_label}.nohydrogen.dcd')

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
        for folder in [self.host_folder, self.na_folder, self.time_folder, self.input_folder,
                       self.charmminp_folder, self.charmmdat_folder, self.ic_folder, self.mat_folder,
                       self.rtficstr_folder, self.datafolder, self.backupfolder, self.scratchfolder]:
            check_dir_exist_and_make(folder)
        
    def initialize_scratch_folders(self):
        for folder in [self.scr_host, self.scr_na, self.scratchfolder]:
            check_dir_exist_and_make(folder)

    def get_refcrd(self, refcrd):
        copyfile(refcrd, self.crd)
        print(f'cp {refcrd} {self.crd}')
        
    def concatenate_trajectory(self, gmx, simu_folder, type_na, start, end):
        na_folder = path.join(simu_folder, self.host, self.type_na)
        roughdir = path.join(na_folder, 'data','roughtrj','1000')
        alltrajfiles = ''
        for mdnum in range(start, end+1):
            filename = path.join(roughdir, '{0}.nopbc.fit.{1}.1000.xtc'.format(type_na, str(mdnum)))
            alltrajfiles = alltrajfiles + filename + ' '
        command = '{0} trjcat -f {1} -o {2} -dt 100'.format(gmx, alltrajfiles, self.xtc_in)
        print(command)
        system(command)

    def remove_redundant_trajectories(self):
        for traj in [self.xtc_in, self.dcd_in]:
            if path.exists(traj):
                info = f'{traj} exists.\n rm {traj}'
                remove(traj)
            else:
                info = f'{traj} dost not exist.'
            print(info)

    def vmd_xtc_to_dcd(self, refpdb):
        vmdscript = path.join(self.input_folder, 'xtc2dcd.vmd')
        self.write_xtc2dcd_vmdscript(vmdscript, refpdb)
        cmd = f'vmd -dispdev text -e {vmdscript}'
        system(cmd)
        print(cmd)

    def write_xtc2dcd_vmdscript(self, vmd_out, refpdb):
        f = open(vmd_out, 'w')
        lines = [f'mol new {refpdb} type pdb',
                 f'mol addfile {self.xtc_in} type xtc waitfor -1',
                 f'animate write dcd {self.dcd_in} beg 1 end 10001 waitfor all',
                 f'exit']
        for line in lines:
            f.write(f'{line}\n')
        f.close()
        print(f'write {vmd_out}')

    def check_size(self, key):
        d_file = {'dcd_in': self.dcd_in, 'dcd_out': self.dcd_out}
        b = path.getsize(d_file[key])
        if b > 0:
            print(f'{self.dcd_in} generate succesfully. Size is {b}.')
        else:
            print(f'{self.dcd_in} Something wrong!!!!! Size is {b}.')

    def make_no_h_dcd_input(self, crd_inp, mass_weighted=True, amber=False, begin=1, frame_num=10000, firstter=None, lastter=None):
        if self.type_na == 'arna+arna':
            supplement1 = None
            supplement2 = None 
        else:
            supplement1 = get_patch(self.seq1, 1)
            supplement2 = get_patch(self.seq2, 2)

        f_inp = path.join(self.charmminp_folder, 'write_no_h_dcd.inp')

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
        inp.read_traj(self.dcd_in)
        inp.open_write_traj(self.dcd_out)
        inp.write_noh_dcd(mass_weighted=mass_weighted, begin=begin, frame_num=frame_num)
        inp.close_unit(21)
        inp.close_unit(30)
        inp.end()

    def make_no_h_dcd(self):
        f_inp = path.join(self.charmminp_folder, 'write_no_h_dcd.inp')
        f_dat = path.join(self.charmmdat_folder, 'write_no_h_dcd.dat')
        exec_charmm(f_inp, f_dat)
