from os import path, system
import numpy as np

class PBSAgent:
    logfolder = '/home/yizaochen/log'

    def __init__(self, filename, jobname, walltime, n_node, n_cpu, pythonexec):
        self.filename = filename
        self.jobname = jobname
        self.walltime = walltime
        self.n_node = n_node
        self.n_cpu = n_cpu
        self.logfile = path.join(self.logfolder, f'{self.jobname}.log')
        self.pythonexec = pythonexec

        self.f = None
        
    def open_file(self):
        self.f = open(self.filename, 'w')

    def close_file(self):
        self.f.close()
        
    def pbssuffix(self):
        self.f.write("#!/bin/bash -l\n")
        self.f.write("#PBS -N {0}\n".format(self.jobname))
        self.f.write("#PBS -l walltime={0}\n".format(self.walltime))
        self.f.write("#PBS -q batch\n")
        self.f.write("#PBS -l nodes={0}:ppn={1}\n".format(self.n_node, self.n_cpu))
        self.f.write("#PBS -j oe\n")
        self.f.write("#PBS -o {0}\n".format(self.logfile))
        self.f.write('#PBS -r n\n')

    def set_pythonexec(self):
        self.f.write('pythonexec=\'{0}\'\n\n'.format(self.pythonexec))
        
    def cutomize_part(self, path_to_py):
        self.f.write(f'program={path_to_py}\n')
        self.f.write('$pythonexec $program\n\n')
        
if __name__ == '__main__':
    host = 'pnas_amber_16mer'
    type_na = 'bdna+bdna'
    start = 0  # default: 0
    end = 200 # default : 400
    cutoff_list = np.arange(6.0, 10.1, 0.5) 
    
    write_qsub = False
    exec_qsub = True
    
    run_n_node = 1
    run_n_cpu = 16
    wtime = '48:00:00'
    qsubscriptfolder = path.join(rootfolder, 'qsub_script')
    
    if write_qsub:
        for cutoff in cutoff_list:
            qsub_fname1 = f'{host}_{type_na}_{cutoff:.2f}.qsub'
            qsub_fname = path.join(qsubscriptfolder, qsub_fname1)
            j_name = '{0}_{1:.2f}'.format(d_abbr[type_na], cutoff)
            logf = path.join(logfolder, '{0}.out'.format(j_name))
            
            p_agent = PBSAgent(qsub_fname, j_name, wtime, run_n_node, run_n_cpu, logf)
            
            p_agent.open_file()
            p_agent.pbssuffix()
            p_agent.rootfolder_exec()
            p_agent.cutomize_part(host, type_na, cutoff, start, end)
            p_agent.close_file()
            
    if exec_qsub:
        for cutoff in cutoff_list:
            qsub_fname1 = f'{host}_{type_na}_{cutoff:.2f}.qsub'
            qsub_fname = path.join(qsubscriptfolder, qsub_fname1)
            cmd = 'qsub {0}'.format(qsub_fname)
            print(cmd)
            system(cmd)    
    
