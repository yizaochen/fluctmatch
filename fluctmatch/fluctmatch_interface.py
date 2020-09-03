from os import path, mkdir
from shutil import copyfile
from subprocess import check_call
from argparse import ArgumentParser
import re
import datetime
import time
from scipy import constants
import numpy as np
from fluctmatch.miscell import check_dir_exist_and_make
from fluctmatch import enm, rtf, ic_str, charmm_fluct, ic_table, prm


T = 310.0 # temperature, 310 K
RT = T * (constants.k * constants.N_A / (constants.calorie * constants.kilo)) # RT kcal/mol # https://en.wikipedia.org/wiki/KT_(energy)
alpha = RT * 0.02  # Learning Rate


def exec_charmm(charmm, f_input, f_output):
    print("charmm< {0} > {1}".format(f_input, f_output))
    check_call(charmm, stdin=open(f_input, 'r'), stdout=open(f_output, 'w+'), shell=True)


def write_ic_fluct_inp(f_out, host, type_na, cutoff, rootfolder):
    inp = charmm_fluct.Script(f_out)
    inp.write_bomlev()
    inp.initialize_variables(rootfolder, host, type_na)
    inp.read_rtf(cutoff)
    inp.read_seq()
    inp.read_crd()
    inp.stream_str(cutoff)
    inp.read_traj()
    inp.icfluct()
    inp.write_icfluct(cutoff)
    inp.end()
    
def write_ic_avg_inp(f_out, host, type_na, cutoff, rootfolder, distance_average=False):
    inp = charmm_fluct.Script(f_out)
    inp.write_bomlev()
    inp.initialize_variables(rootfolder, host, type_na)
    inp.read_rtf(cutoff)
    inp.read_seq()
    inp.read_crd()
    inp.stream_str(cutoff)
    inp.read_traj()
    inp.icavg(distance_average=distance_average)
    inp.write_icavg(cutoff)
    inp.end()
    
def write_nmainit_inp(f_out, host, type_na, cutoff, rootfolder, scratchfolder):
    inp = charmm_fluct.Script(f_out)
    inp.write_bomlev()
    inp.initialize_variables_nma(rootfolder, scratchfolder, host, type_na, cutoff)
    inp.read_rtf(cutoff)
    inp.read_prm()
    inp.read_seq()
    inp.read_crd_nma()
    inp.set_mass_1()
    inp.minimization()
    inp.stream_str_nma()
    inp.nma()
    inp.end()
    
    
def fluct_match(host, type_na, start, end, cutoff, icfluct_0, icavg_0, scratchfolder, charmm, rootfolder):
    # Backup Folder
    backupfolder = path.join(scratchfolder, 'backup')
    check_dir_exist_and_make(backupfolder)
    
    # Set initial variables
    avg_file = path.join(scratchfolder, 'average_{0:.2f}.ic'.format(cutoff))
    fluct_file = path.join(scratchfolder, 'fluct_{0:.2f}.ic'.format(cutoff))
    prmfile = path.join(scratchfolder, 'na_enm_{0:.2f}.prm'.format(cutoff)) 
    prm_backup = path.join(backupfolder, 'na_enm_{0:.2f}.backup.prm'.format(cutoff))
    prm_backup_0 = path.join(backupfolder, 'na_enm_{0:.2f}.backup.0.prm'.format(cutoff))
    err_file = path.join(scratchfolder, 'error_{0:.2f}.txt'.format(cutoff))
    copyfile(prmfile, prm_backup_0)
    
    # Read icavg, icfluct, set k and target
    icavg = ic_table.ICTable(avg_file)
    icfluct = ic_table.ICTable(fluct_file)
    kbpair = ic_table.KBPair(read_from_prm=False, icavg=icavg, icfluct=icfluct, rt=RT)
    k = kbpair.d['k']
    target = np.reciprocal(np.square(icfluct_0.values))
    
    # Write Error
    f = open(err_file, 'w')
    f.write('Created at {0}\n'.format(datetime.datetime.now()))
    f.write('{0:<5} {1:8}\n'.format('n_iter', 'error'))
    f.close()
   
    # Iteration and Do NMA
    for i in range(start, end):
        print("IterNum: {0}".format(i))
        nma_inp = path.join(scratchfolder, 'nma.inp')
        nma_dat = path.join(scratchfolder, 'nma.dat')
        if i == 0:
            write_nmainit_inp(nma_inp, host, type_na, cutoff, rootfolder, scratchfolder)
        exec_charmm(charmm, nma_inp, nma_dat)
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
        
        # Update PRM
        prm_agent = prm.PRM(host, type_na, new_kbpair, iternum=i)
        prm_agent.write_prm(prmfile)


def main(rootfolder, host, type_na, cutoff, start, end, charmm):
    # Initialize 
    agent = enm.ENMAgent(rootfolder, host, type_na)
    nadir = agent.na_folder
    charmminpfolder = path.join(nadir, 'charmm_inp')
    charmmdatfolder = path.join(nadir, 'charmm_dat')
    icfolder = path.join(nadir, 'ic')
    datafolder = path.join(nadir, 'data')
    scratchfolder = path.join(nadir, 'scratch') 

    # IC fluct
    icfluct_inp = path.join(charmminpfolder, 'ic_fluct_{0:.2f}.inp'.format(cutoff))
    icfluct_dat = path.join(charmmdatfolder, 'ic_fluct_{0:.2f}.dat'.format(cutoff))
    write_ic_fluct_inp(icfluct_inp, host, type_na, cutoff, rootfolder)
    exec_charmm(charmm, icfluct_inp, icfluct_dat)
    mode0ic = path.join(icfolder, 'mode.0.{0:.2f}.ic'.format(cutoff))
    nafluctic = path.join(datafolder, 'na.fluct.{0:.2f}.ic'.format(cutoff))
    with open(mode0ic, 'r') as f:
        context = f.read()
    context = re.sub(r'-99 ', ' -99 ', context)
    with open(nafluctic, 'w') as f:
        f.write(context)
    icfluct_0 = ic_table.ICTable(nafluctic, initial=True)

    # IC Avg
    icavg_inp = path.join(charmminpfolder, 'ic_avg_{0:.2f}.inp'.format(cutoff))
    icavg_dat = path.join(charmmdatfolder, 'ic_avg_{0:.2f}.dat'.format(cutoff))
    write_ic_avg_inp(icavg_inp, host, type_na, cutoff, rootfolder, distance_average=False) # Important! Check Fix b0
    exec_charmm(charmm, icavg_inp, icavg_dat)
    mode0avgic = path.join(icfolder, 'mode.0.avg.{0:.2f}.ic'.format(cutoff))
    naavgic = path.join(datafolder, 'na.avg.{0:.2f}.ic'.format(cutoff))
    with open(mode0avgic, 'r') as f:
        context = f.read()
    context = re.sub(r'-99 ', ' -99 ', context)
    with open(naavgic, 'w') as f:
        f.write(context)
    icavg_0 = ic_table.ICTable(naavgic, initial=True)

    # Get the initial equilibrium distance and force constant and write PRM
    b_0 = icavg_0.values  # Initial Guess of equilibrium bond length
    k_0 = RT / np.square(icfluct_0.values) # Initial Guess of force constants
    kbpair_0 = ic_table.KBPair(read_from_prm=False, icavg=icavg_0, icfluct=icfluct_0, rt=RT)
    scratch_prm = path.join(scratchfolder, 'na_enm_{0:.2f}.prm'.format(cutoff))  #!!! Important File
    prm_agent = prm.PRM(host, type_na, kbpair_0, iternum=0)
    prm_agent.write_prm(scratch_prm)

    # NMA Initialize
    nma_init_inp = path.join(scratchfolder, 'nmainit.inp')
    nma_init_dat = path.join(scratchfolder, 'nmainit.dat')
    write_nmainit_inp(nma_init_inp, host, type_na, cutoff, rootfolder, scratchfolder)
    exec_charmm(charmm, nma_init_inp, nma_init_dat)
    
    # Fluctuation-Matching
    fluct_match(host, type_na, start, end, cutoff, icfluct_0, icavg_0, scratchfolder, charmm, rootfolder)
    
    # Copy Back to the folder
    cutoffdatafolder = path.join(nadir, 'cutoffdata')
    check_dir_exist_and_make(cutoffdatafolder)
    lastprmfile = path.join(cutoffdatafolder, 'na_enm_{0:.2f}.prm'.format(cutoff))
    copyfile(scratch_prm, lastprmfile)
    print('cp {0} {1}'.format(scratch_prm, lastprmfile))
    
    err_file = path.join(scratchfolder, 'error_{0:.2f}.txt'.format(cutoff))
    lasterrfile = path.join(cutoffdatafolder, 'error_{0:.2f}.txt'.format(cutoff))
    copyfile(err_file, lasterrfile)
    print('cp {0} {1}'.format(err_file, lasterrfile))
    
    nma_inp = path.join(scratchfolder, 'nma.inp')
    nma_inp_back = path.join(nadir, 'charmm_inp', f'nma_{cutoff:.2f}.inp')
    copyfile(nma_inp, nma_inp_back)
    
    nma_dat = path.join(scratchfolder, 'nma.dat')
    nma_dat_back = path.join(nadir, 'charmm_dat', f'nma_{cutoff:.2f}.dat')
    copyfile(nma_dat, nma_dat_back)

