from os import path, mkdir
from shutil import copyfile
from subprocess import check_call
import re
import datetime
import time
from scipy import constants
import numpy as np
import enm, rtf, ic_str, charmm_fluct, ic_table, prm

charmm = '/Users/yizao/c41b1_yz/exec/osx/charmm'
scratchfolder = '/scratch/yizaochen' 

T = 310.0 # temperature, 310 K
RT = T * (constants.k * constants.N_A / (constants.calorie * constants.kilo)) # RT kcal/mol # https://en.wikipedia.org/wiki/KT_(energy)
alpha = RT * 0.02  # Learning Rate

def check_dir_exist_and_make(file_path):
    if path.exists(file_path):
        print("{0} exists".format(file_path))
    else:
        print("mkdir {0}".format(file_path))
        mkdir(file_path)


def exec_charmm(f_input, f_output):
    print("charmm< {0} > {1}".format(f_input, f_output))
    check_call(charmm, stdin=open(f_input, 'r'), stdout=open(f_output, 'w+'), shell=True)

    
"""
The process of Fluctuation-Matching
"""
def write_ic_fluct_inp(f_out, host, type_na):
    inp = charmm_fluct.Script(f_out)
    inp.write_bomlev()
    inp.initialize_variables(host, type_na)
    inp.read_rtf()
    inp.read_seq()
    inp.read_crd()
    inp.stream_str()
    inp.read_traj()
    inp.icfluct()
    inp.write_icfluct()
    inp.end()
    
    
def write_ic_avg_inp(f_out, host, type_na, distance_average=False):
    inp = charmm_fluct.Script(f_out)
    inp.write_bomlev()
    inp.initialize_variables(host, type_na)
    inp.read_rtf()
    inp.read_seq()
    inp.read_crd()
    inp.stream_str()
    inp.read_traj()
    inp.icavg(distance_average=distance_average)
    inp.write_icavg()
    inp.end()
    
    
def write_nmainit_inp(f_out, host, type_na, out_start_end_mode=None):
    inp = charmm_fluct.Script(f_out)
    inp.write_bomlev()
    inp.initialize_variables_nma(host, type_na)
    inp.read_rtf()
    inp.read_prm()
    inp.read_seq()
    inp.read_crd_nma()
    inp.set_mass_1()
    inp.minimization()
    inp.stream_str_nma()
    inp.nma(out_start_end_mode=out_start_end_mode)
    inp.end()
    
    
def write_nma_inp(f_out, host, type_na, out_start_end_mode=None):
    inp = charmm_fluct.Script(f_out)
    inp.write_bomlev()
    inp.initialize_variables_nma(host, type_na)
    inp.read_rtf()
    inp.read_prm()
    inp.read_seq()
    inp.read_crd_nma()
    inp.minimization()
    inp.stream_str_nma()
    inp.nma(out_start_end_mode=out_start_end_mode)
    inp.end()
    
    
def fluct_match(host, type_na, start, end, icfluct_0, icavg_0, kbpair_0, nadir, out_start_end_mode=None):
    # Important folders
    charmminpfolder = path.join(nadir, 'charmm_inp')
    charmmdatfolder = path.join(nadir, 'charmm_dat')
    datafolder = path.join(nadir, 'data')
    iterfolder = path.join(nadir, 'diff_iters')
    check_dir_exist_and_make(iterfolder)
    backupfolder = path.join(datafolder, 'backup')
    check_dir_exist_and_make(backupfolder)
    
    # Set initial variables
    avg_file = path.join(datafolder, 'average.ic')
    fluct_file = path.join(datafolder, 'fluct.ic')
    prmfile = path.join(datafolder, 'na_enm.prm') 
    prm_backup = path.join(backupfolder, 'na_enm.backup.prm')
    prm_backup_0 = path.join(backupfolder, 'na_enm.backup.0.prm')
    err_file = path.join(datafolder, 'error.txt')
    copyfile(prmfile, prm_backup_0)
    
    # Read icavg, icfluct, set k and target
    icavg = ic_table.ICTable(avg_file)
    icfluct = ic_table.ICTable(fluct_file)
    kbpair = ic_table.KBPair(read_from_prm=False, icavg=icavg, icfluct=icfluct, rt=RT)
    k = kbpair_0.d['k']
    target = np.reciprocal(np.square(icfluct_0.values))
    
    # Write Error
    f = open(err_file, 'w')
    f.write('Created at {0}\n'.format(datetime.datetime.now()))
    f.write('{0:<5} {1:8}\n'.format('n_iter', 'error'))
    f.close()
   
    # Iteration and Do NMA
    for i in range(start, end+1):
        print("IterNum: {0}".format(i))
        nma_inp = path.join(charmminpfolder, 'nma.inp')
        nma_dat = path.join(charmmdatfolder, 'nma.dat')
        if i == 0:
            write_nmainit_inp(nma_inp, host, type_na, out_start_end_mode=out_start_end_mode)
        exec_charmm(nma_inp, nma_dat)
        time.sleep(2)
        
        # Read icavg, icfluct
        new_icavg = ic_table.ICTable(avg_file)
        new_icfluct = ic_table.ICTable(fluct_file)
        optimized = np.reciprocal(np.square(new_icfluct.values))
        # Be Careful, b0 is from icavg_0, which is the average structure of MD, it should be kept fix
        new_kbpair = ic_table.KBPair(read_from_prm=False, icavg=icavg_0, icfluct=new_icfluct, rt=RT)

        # Update k
        k -= alpha * (optimized - target)
        k[np.where(k < 0.)] = 0. # Make all negative values zero
        new_kbpair.d['k'] = k
        
        # Calculate RMSD error and Write
        error = icfluct_0.values - new_icfluct.values
        error = np.square(error).mean()
        error = np.sqrt(error)
        f = open(err_file, 'a')
        f.write('{0:<5} {1:8.4f}\n'.format(i+1, error))
        f.close()
        
        # Backup PRM
        copyfile(prmfile, prm_backup)
        
        # Copy prm with different iteration number to original folder
        iternum_prm = path.join(iterfolder, f'na_enm_{i}.prm')
        copyfile(prmfile, iternum_prm)
        
        # Update PRM
        prm_agent = prm.PRM(host, type_na, new_kbpair, iternum=i)
        prm_agent.write_prm(prmfile)
        
        
"""
Post-Process after Fluctuation-Matching
"""
def write_mini_inp(f_out, prm_in, host, type_na, cutoff, crd_out):
    inp = charmm_fluct.Script(f_out)
    inp.initialize_variables_mini(host, type_na, prm_in)
    inp.read_rtf_cutoff(cutoff)
    inp.read_prm()
    inp.read_seq()
    inp.read_crd_nma()
    inp.set_mass_1()
    inp.minimization_mini()
    inp.write_crd_mini_by_filename(crd_out)
    inp.end()
    print(f'make charmm inp file: {f_out}')
        
        
"""
Functions for universal k
"""
def make_universal_k_prm_based_on_fm_result(host, type_na, cutoff, clean_criteria, k_value):
    """
    Read a prm which is a result of fluctuation-matching, only make those remained force constants a universal value
    """
    prm_in = path.join('/Users/yizao/PycharmProjects/ENM', host, type_na, 'diffcutoff', f'na_enm_{cutoff:.2f}.prm')
    kbpair = ic_table.KBPair(read_from_prm=True, filename=prm_in)
    closezeros_indices = np.argwhere(kbpair.d['k'] < clean_criteria)
    closezeros_indices = np.stack(closezeros_indices, axis=1)[0]
    kbpair.d['k'][closezeros_indices] = 0
    nonzero_indices = np.nonzero(kbpair.d['k'])[0]
    n_nonzero_bonds = len(nonzero_indices)
    print(f'There are {n_nonzero_bonds} non-zero k bonds.')
    kbpair.d['k'][nonzero_indices] = k_value
    prm_out = path.join('/Users/yizao/PycharmProjects/ENM', host, type_na, 'diffcutoff', f'na_enm_{cutoff:.2f}_universalk_afterfm.prm')
    prm_agent = prm.PRM(host, type_na, kbpair, iternum=0)
    prm_agent.write_prm(prm_out)
    print(f'Write PRM file: {prm_out}')
    return prm_out


def make_universal_k_prm_all_initialbonds(host, type_na, cutoff, clean_criteria, k_value):
    """
    Let all bonds are universal value, these bonds are created when building ENM initially
    """
    prm_in = path.join('/Users/yizao/PycharmProjects/ENM', host, type_na, 'diffcutoff', f'na_enm_{cutoff:.2f}.prm')
    kbpair = ic_table.KBPair(read_from_prm=True, filename=prm_in)
    kbpair.d['k'][:] = k_value
    nonzero_indices = np.nonzero(kbpair.d['k'])[0]
    n_nonzero_bonds = len(nonzero_indices)
    print(f'There are {n_nonzero_bonds} non-zero k bonds.')
    prm_out = path.join('/Users/yizao/PycharmProjects/ENM', host, type_na, 'diffcutoff', f'na_enm_{cutoff:.2f}_universalk_initial_enmbonds.prm')
    prm_agent = prm.PRM(host, type_na, kbpair, iternum=0)
    prm_agent.write_prm(prm_out)
    print(f'Write PRM file: {prm_out}')
    return prm_out
    