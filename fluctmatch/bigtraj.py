from os import path, system, remove
from shutil import copyfile, make_archive
import zipfile
from fluctmatch.avg_dcd_noh import AvgcrddcdAgent, exec_charmm
from fluctmatch.miscell import check_dir_exist_and_make, get_patch
from fluctmatch.sequence import sequences
from fluctmatch.charmm import Script
from fluctmatch.enm import ENMAgent
from fluctmatch.rtf import RTF
from fluctmatch.ic_str import ICSTR
from fluctmatch.fluct_pbs import PBSAgent

class BigTrajAgent(AvgcrddcdAgent):

    start_time = 0
    end_time = 5000 # 5000 ns
    interval_time = 1000
    gmx = '/usr/bin/gmx'
    
    def __init__(self, host, type_na, allsys_folder, bigtraj_folder, simu_folder, split_5=True, one_big_window=False):
        super().__init__(host, type_na, allsys_folder)
        self.allsys_folder = allsys_folder
        self.bigtraj_folder = bigtraj_folder
        self.simu_folder = simu_folder

        if split_5:
            self.time_list, self.mdnum_list = self.get_time_list_split_5()
            self.multiscale_bigtraj_folder = '/home/yizaochen/bigtraj_fluctmatch/split_5'
        elif one_big_window:
            self.time_list, self.mdnum_list = self.get_time_list_one_big_window()
            self.multiscale_bigtraj_folder = '/home/yizaochen/bigtraj_fluctmatch/5000ns'
        else:
            self.time_list, self.mdnum_list = self.get_time_list()
            self.multiscale_bigtraj_folder = '/home/yizaochen/bigtraj_fluctmatch/1000ns'

        self.d_smallagents = self.get_all_small_agents()

        self.d_pairs = None

    def get_time_list_split_5(self):
        n_split = 5 # ad hoc
        mdnum_list = list()
        time_list = list()
        for time1 in range(n_split):
            time2 = time1 + 1
            time_list.append((time1, time2))
        return time_list, mdnum_list   

    def get_time_list_one_big_window(self):
        mdnum_list = list()
        time_list = [(0, 5000)]
        return time_list, mdnum_list     

    def get_time_list(self):
        middle_interval = int(self.interval_time/2)
        time_list = list()
        mdnum_list = list()
        mdnum1 = 1
        for time1 in range(self.start_time, self.end_time, middle_interval):
            time2 = time1 + self.interval_time
            if time2 <= self.end_time:
                time_list.append((time1, time2))
                mdnum_list.append((mdnum1, mdnum1+9)) # This should be more general
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

    def copy_no_h_5us_dcd_to_smallfolder(self):
        for time1, time2 in self.time_list:
            noh_dcd = path.join(self.allsys_folder, self.host, self.type_na, 'input', 'heavyatoms', f'{self.type_na}.nohydrogen.dcd')
            self.d_smallagents[(time1,time2)].get_dcdout(noh_dcd)

    def check_no_h_5us_dcd(self):
        for time1, time2 in self.time_list:
            self.d_smallagents[(time1,time2)].vmd_check_dcdout()

    def concatenate_xtc_by_gmx_split_5(self, mdnum1=1, mdnum2=50):
        start_time = 100
        for time1, time2 in self.time_list:
            self.d_smallagents[(time1,time2)].concatenate_trajectory_split_5(self.gmx, self.simu_folder, self.type_na, mdnum1, mdnum2, start_time)
            start_time += 100

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

    def make_avg_crd(self):
        for time1, time2 in self.time_list:
            self.d_smallagents[(time1,time2)].make_avg_crd_input(amber=True, firstter='amber_5ter', lastter='amber_3ter')
            self.d_smallagents[(time1,time2)].make_avg_crd()

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

    def update_to_multiphysics(self):
        old_f = path.join(self.bigtraj_folder, f'{self.host}.zip')
        new_f = path.join(self.multiscale_bigtraj_folder, f'{self.host}.zip')
        cmd = f'scp {old_f} yizaochen@multiphysics:{new_f}'
        print(cmd)

    def download_resultzip_from_multiphysics(self):
        old_f = path.join(self.multiscale_bigtraj_folder, 'zipfiles', f'{self.host}_results.zip')
        new_f = path.join(self.bigtraj_folder, f'{self.host}_results.zip')
        cmd = f'scp yizaochen@multiphysics:{old_f} {new_f}'
        print(cmd)

    def download_resultzip_from_multiscale(self):
        old_f = path.join(self.multiscale_bigtraj_folder, 'zipfiles', f'{self.host}_results.zip')
        new_f = path.join(self.bigtraj_folder, f'{self.host}_results.zip')
        cmd = f'scp yizaochen@multiscale:{old_f} {new_f}'
        print(cmd)

    def unzip_to_tempresults(self):
        temp_results = path.join(self.bigtraj_folder, 'temp_results')
        temp_host_folder = path.join(temp_results, self.host)
        f_zip = path.join(self.bigtraj_folder, f'{self.host}_results.zip')
        check_dir_exist_and_make(temp_host_folder)
        with zipfile.ZipFile(f_zip, 'r') as zip_ref:
            zip_ref.extractall(temp_host_folder)
        info = f'Unzip {f_zip} into {temp_host_folder}'
        print(info)

    def redistribute_results(self, cutoff):
        temp_results = path.join(self.bigtraj_folder, 'temp_results')
        temp_host_folder = path.join(temp_results, self.host)
        for time1, time2 in self.time_list:
            agent = self.d_smallagents[(time1,time2)]
            agent.redistribute(temp_host_folder, cutoff)      

class BigTrajOnServer(BigTrajAgent):
    def __init__(self, host, type_na, bigtraj_folder, split_5=True, one_big_window=False):
        self.host = host
        self.type_na = type_na
        self.bigtraj_folder = bigtraj_folder

        self.host_folder = path.join(bigtraj_folder, host)
        self.na_folder = path.join(self.host_folder, type_na)
        
        self.zipfolder = path.join(self.bigtraj_folder, 'zipfiles')
        self.f_zip = path.join(self.zipfolder, f'{host}.zip')

        if split_5:
            self.time_list, self.mdnum_list = self.get_time_list_split_5()
        elif one_big_window:
            self.time_list, self.mdnum_list = self.get_time_list_one_big_window()
        else:
            self.time_list, self.mdnum_list = self.get_time_list()

        self.d_smallagents = self.get_all_small_agents()

    def unarchive_folder(self):
        check_dir_exist_and_make(self.host_folder)
        with zipfile.ZipFile(self.f_zip, 'r') as zip_ref:
            zip_ref.extractall(self.host_folder)
        info = f'Unzip {self.f_zip} into {self.host_folder}'
        print(info)

    def generate_python_files(self, start, end):
         for time1, time2 in self.time_list:
            agent = self.d_smallagents[(time1,time2)]
            agent.make_python_fluct_main(start, end)

    def generate_qsub_scripts(self, path_to_pythonexec):
        for time1, time2 in self.time_list:
            agent = self.d_smallagents[(time1,time2)]
            agent.make_qsub_file(path_to_pythonexec)

    def submit_all_qsubs(self):
        for time1, time2 in self.time_list:
            agent = self.d_smallagents[(time1,time2)]
            agent.submit_qsub()

    def zip_all_results(self, cutoff):
        for time1, time2 in self.time_list:
            agent = self.d_smallagents[(time1,time2)]
            agent.check_result_zipfolder(self.zipfolder)
            agent.copy_result_to_zipfolder(self.zipfolder, cutoff)
        output_name = path.join(self.zipfolder, f'{self.host}_results')
        target_name = path.join(self.zipfolder, self.host)
        make_archive(output_name, 'zip', target_name)
        print(f'Archive {target_name} into {output_name}.zip')


class SmallTrajAgent(ENMAgent):
    scratchroot = '/scratch'

    def __init__(self, root, host, type_na, time_label):
        self.rootfolder = root
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
        self.qsub_folder = path.join(self.time_folder, 'qsub_scripts')
        self.qsub_file = path.join(self.qsub_folder, 'fluct_main.qsub')
        self.pyfile = path.join(self.time_folder, 'fluct_main.py')
        #self.initialize_folders()

        # Scratch folder
        self.scr_host = path.join(self.scratchroot, host)
        self.scr_na = path.join(self.scr_host, type_na)
        self.scratchfolder = path.join(self.scr_na, self.time_label)

        self.seq1 = sequences[self.host][self.type_na]['guide']
        self.seq2 = sequences[self.host][self.type_na]['target']

        self.enmprm = path.join(self.datafolder, 'na_enm.prm')
        #self.mode0dcd = path.join(self.mode_traj_folder, 'mode.0.dcd')
        self.avg_crd = path.join(self.input_folder, '{0}.nohydrogen.avg.crd'.format(self.type_na))
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
                       self.rtficstr_folder, self.datafolder, self.backupfolder]:
            check_dir_exist_and_make(folder)
        
    def initialize_scratch_folders(self):
        for folder in [self.scr_host, self.scr_na, self.scratchfolder]:
            check_dir_exist_and_make(folder)

    def get_refcrd(self, refcrd):
        copyfile(refcrd, self.crd)
        print(f'cp {refcrd} {self.crd}')

    def get_dcdout(self, noh_dcd):
        copyfile(noh_dcd, self.dcd_out)
        print(f'cp {noh_dcd} {self.dcd_out}')

    def vmd_check_dcdout(self):
        print(f'vmd -cor {self.crd} {self.dcd_out}')

    def concatenate_trajectory_split_5(self, gmx, simu_folder, type_na, start, end, starttime):
        na_folder = path.join(simu_folder, self.host, self.type_na)
        roughdir = path.join(na_folder, 'data','roughtrj','1000')
        alltrajfiles = ''
        for mdnum in range(start, end+1):
            filename = path.join(roughdir, '{0}.nopbc.fit.{1}.1000.xtc'.format(type_na, str(mdnum)))
            alltrajfiles = alltrajfiles + filename + ' '
        command = '{0} trjcat -f {1} -o {2} -dt 100 -b {3}'.format(gmx, alltrajfiles, self.xtc_in, starttime)
        print(command)
        system(command)
        
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

    def make_avg_crd_input(self, selection='all', amber=False, firstter=None, lastter=None):
        if self.type_na == 'arna+arna':
            supplement1 = None
            supplement2 = None 
        else:
            supplement1 = get_patch(self.seq1, 1)
            supplement2 = get_patch(self.seq2, 2)

        f_inp = path.join(self.charmminp_folder, 'write_no_h_avg_crd.inp')

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
        inp.read_crd(self.crd)
        inp.read_traj(self.dcd_out)
        inp.calculate_avg(selection=selection)
        inp.coor_copy()
        inp.write_crd(self.avg_crd, comp=True)
        inp.end()

    def make_avg_crd(self):
        f_inp = path.join(self.charmminp_folder, 'write_no_h_avg_crd.inp')
        f_dat = path.join(self.charmmdat_folder, 'write_no_h_avg_crd.dat')
        exec_charmm(f_inp, f_dat)

    def make_python_fluct_main(self, start, end):
        lines = ['from fluctmatch import fluctmatch_interface\n',
                 f'bigtraj_folder = \'{self.rootfolder}\'',
                 f'charmm = \'/home/yizaochen/c39b1_yizao/exec/gnu/charmm\'',
                 f'host = \'{self.host}\'',
                 f'type_na = \'{self.type_na}\'',
                 f'time_label = \'{self.time_label}\'',
                 'cutoff=4.7',
                 f'start={start}',
                 f'end={end}\n',
                 f'fluctmatch_interface.main_split_window(bigtraj_folder, host, type_na, time_label, cutoff, start, end, charmm)']
        f = open(self.pyfile, 'w')
        for line in lines:
            f.write(line)
            f.write('\n')
        f.close()
        print(f'Generate {self.pyfile}')

    def make_qsub_file(self, path_to_pythonexec):
        check_dir_exist_and_make(self.qsub_folder)
        jobname = f'{self.host}_{self.time_label}'
        wtime = '48:00:00'
        run_n_node = 1
        run_n_cpu = 16
        p_agent = PBSAgent(self.qsub_file, jobname, wtime, run_n_node, run_n_cpu, path_to_pythonexec)
        p_agent.open_file()
        p_agent.pbssuffix()
        p_agent.set_pythonexec()
        p_agent.cutomize_part(self.pyfile)
        p_agent.close_file()

    def submit_qsub(self):
        cmd = 'qsub {0}'.format(self.qsub_file)
        print(cmd)
        system(cmd)  

    def check_result_zipfolder(self, ziproot):
        zip_host_folder = path.join(ziproot, self.host)
        zip_na_folder = path.join(zip_host_folder, self.type_na)
        zip_time_folder = path.join(zip_na_folder, self.time_label)
        for folder in [zip_host_folder, zip_na_folder, zip_time_folder]:
            check_dir_exist_and_make(folder)

    def copy_result_to_zipfolder(self, ziproot, cutoff):
        zip_host_folder = path.join(ziproot, self.host)
        zip_na_folder = path.join(zip_host_folder, self.type_na)
        zip_time_folder = path.join(zip_na_folder, self.time_label)
        oldfile = path.join(self.time_folder, 'cutoffdata', f'na_enm_{cutoff:.2f}.prm')
        newfile = path.join(zip_time_folder, f'na_enm_{cutoff:.2f}.prm')
        copyfile(oldfile, newfile)
        print(f'cp {oldfile} {newfile}')

    def redistribute(self, temp_host_folder, cutoff):
        oldfile = path.join(temp_host_folder, self.type_na, self.time_label, f'na_enm_{cutoff:.2f}.prm')
        newfile = path.join(self.datafolder, f'na_enm_{cutoff:.2f}.prm')
        copyfile(oldfile, newfile)
        print(f'cp {oldfile} {newfile}')       