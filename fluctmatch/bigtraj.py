from os import path
from fluctmatch.avg_dcd_noh import AvgcrddcdAgent
from fluctmatch.miscell import check_dir_exist_and_make

class BigTrajAgent(AvgcrddcdAgent):

    start_time = 0
    end_time = 5000 # 5000 ns
    interval_time = 1000 # unit: ns

    def __init__(self, host, type_na, allsys_folder, bigtraj_folder, simu_folder):
        super().__init__(host, type_na, allsys_folder)
        self.allsys_folder = allsys_folder
        self.bigtraj_folder = bigtraj_folder
        self.simu_folder = simu_folder

        self.time_list = self.get_time_list()

    def get_time_list(self):
        middle_interval = int(self.interval_time/2)
        time_list = list()
        for time1 in range(self.start_time, self.end_time, middle_interval):
            time2 = time1 + self.interval_time
            if time2 <= self.end_time:
                time_list.append((time1, time2))
        return time_list

class SmallTrajAgent:
    def __init__(self, root, host, type_na, time_label):
        self.host = host
        self.type_na = type_na
        self.time_label = time_label

        self.host_folder = path.join(root, host)
        self.na_folder = path.join(self.host_folder, type_na, time_label)
        self.input_folder = path.join(self.na_folder, 'input')
        self.charmminp_folder = path.join(self.na_folder, 'charmm_inp')
        self.charmmdat_folder = path.join(self.na_folder, 'charmm_dat')
        self.mode_traj_folder = path.join(self.na_folder, 'mode_traj')
        self.ic_folder = path.join(self.na_folder, 'ic')
        self.mat_folder = path.join(self.na_folder, 'ic_fluct_mat')
        self.rtficstr_folder = path.join(self.na_folder, 'rtf_ic_str')
        self.datafolder = path.join(self.na_folder, 'data')
        self.backupfolder = path.join(self.datafolder, 'backup')
        self.scratchfolder = path.join(self.na_folder, 'scratch')
        self.initialize_folders()

    def initialize_folders(self):
        for folder in [self.host_folder, self.na_folder, self.input_folder,
                       self.charmminp_folder, self.charmmdat_folder,
                       self.mode_traj_folder, self.ic_folder, self.mat_folder,
                       self.rtficstr_folder, self.datafolder, self.backupfolder, self.scratchfolder]:
            check_dir_exist_and_make(folder)