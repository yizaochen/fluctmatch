from os import path
import numpy as np
import charmm_fluct
from fluctmatch_interface import exec_charmm

rootfolder = '/home/yizaochen/fluct_diffcutoff'

def write_mini_inp(f_out, host, type_na, cutoff):
    inp = charmm_fluct.Script(f_out)
    inp.initialize_variables_mini(host, type_na, cutoff)
    inp.read_rtf(cutoff)
    inp.read_prm()
    inp.read_seq()
    inp.read_crd_nma()
    inp.set_mass_1()
    inp.minimization_mini()
    inp.write_crd_mini(cutoff)
    inp.end()

if __name__ == '__main__':
    host = 'pnas_amber_16mer'
    type_na = 'arna+arna'
    cutoff_list = np.arange(6.0, 10.1, 0.5) 
    
    nadir = path.join(rootfolder, host, type_na)
    charmminpfolder = path.join(nadir, 'charmm_inp')
    charmmdatfolder = path.join(nadir, 'charmm_dat')
    
    for cutoff in cutoff_list:
        mini_inp = path.join(charmminpfolder, f'mini_{cutoff:.2f}.inp')
        mini_dat = path.join(charmmdatfolder, f'mini_{cutoff:.2f}.dat')
        write_mini_inp(mini_inp, host, type_na, cutoff)
        exec_charmm(mini_inp, mini_dat)
        
    