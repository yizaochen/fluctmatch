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
        
    def initialize_variables(self, rootfolder, host, type_na):
        self.f.write(f'set rootfolder {rootfolder} \n')
        self.f.write('set host {0} \n'.format(host))
        self.f.write('set typena {0} \n'.format(type_na))
        self.f.write('set nadir @rootfolder/@host/@typena \n')
        self.f.write('set datadir @nadir/data \n')
        self.f.write('set dcddir  @nadir/mode_traj \n')
        self.f.write('set topdir @nadir/rtf_ic_str \n')
        self.f.write('set modeid = 0 \n\n')
        
    def initialize_variables_nma(self, rootfolder, scratchfolder, host, type_na, cutoff):
        self.f.write(f'set rootfolder {rootfolder} \n')
        self.f.write(f'set scratchfolder {scratchfolder} \n')
        self.f.write('set host {0} \n'.format(host))
        self.f.write('set typena {0} \n'.format(type_na))
        self.f.write('set nadir @rootfolder/@host/@typena \n')
        self.f.write('set datadir @nadir/data \n')
        self.f.write('set dcddir  @nadir/mode_traj \n')
        self.f.write('set topdir @nadir/rtf_ic_str \n')
        self.f.write('set icstr   @topdir/na_enm_{0:.2f}.str \n'.format(cutoff))
        self.f.write('set psffile @nadir/input/na_enm.psf \n')
        self.f.write('set crdfile @nadir/input/na_enm.crd \n')
        self.f.write('set prmfile @scratchfolder/na_enm_{0:.2f}.prm \n'.format(cutoff))
        self.f.write('set vibcrd  @scratchfolder/na_enm.vib_{0:.2f}.crd \n'.format(cutoff))
        self.f.write('set avgic   @scratchfolder/average_{0:.2f}.ic \n'.format(cutoff))
        self.f.write('set fluctic @scratchfolder/fluct_{0:.2f}.ic \n'.format(cutoff))
        self.f.write('set vibic   @scratchfolder/na_enm_{0:.2f}.vib \n\n'.format(cutoff))
        self.f.write('set fileu  10 \n')
        self.f.write('set temp   310.0 \n\n')
        
    def initialize_variables_mini(self, rootfolder, host, type_na, cutoff):
        self.f.write(f'set rootfolder {rootfolder} \n')
        self.f.write('set host {0} \n'.format(host))
        self.f.write('set typena {0} \n'.format(type_na))
        self.f.write('set nadir @rootfolder/@host/@typena \n')
        self.f.write('set datadir @nadir/data \n')
        self.f.write('set topdir @nadir/rtf_ic_str \n')
        self.f.write(f'set prmfile @nadir/cutoffdata/na_enm_{cutoff:.2f}.prm \n')
        self.f.write('set crdfile @nadir/input/na_enm.crd \n')
        self.f.write('set fileu  10 \n\n')
        
    def read_rtf(self, cutoff):
        self.f.write('open unit 11 read form name @topdir/na_enm_{0:.2f}.rtf \n'.format(cutoff))
        self.f.write('read rtf card unit 11 \n')
        self.f.write('close unit 11 \n\n')
        
    def read_prm(self):
        self.f.write('open read unit @fileu form name @prmfile \n')
        self.f.write('read para unit @fileu card \n')
        self.f.write('close unit @fileu \n\n')
        
    def read_seq(self):   
        self.f.write('read sequ card \n')
        self.f.write('* \n')
        self.f.write(' 1\n')
        self.f.write('NA\n')
        self.f.write('generate DNA setup\n\n')
        
    def read_crd(self):
        self.f.write('open read unit 12 card name @nadir/input/na_enm.crd \n')
        self.f.write('read coor unit 12 ignore \n')
        self.f.write('close unit 12 \n')
        
    def read_crd_nma(self):
        self.f.write('open read unit @fileu card name @crdfile \n')
        self.f.write('read coor unit @fileu ignore \n')
        self.f.write('coor copy comp \n')
        self.f.write('close unit @fileu \n\n')
        
    def set_mass_1(self):
        self.f.write('scalar mass set 1 sele all end\n\n')
        
    def minimization(self):
        self.f.write('skip all excl bond \n')
        self.f.write('update inbfrq 0 1hbfrq 0 imgfrq 0 \n')
        self.f.write('ener \n')
        self.f.write('mini sd nstep 100 \n')
        self.f.write('mini abnr nstep 2000 \n\n')
        self.f.write('coor orie rms mass \n')
        self.f.write('scalar wmain copy mass \n\n')
        self.f.write('ioformat extended \n')
        self.f.write('open write unit @fileu card name @vibcrd \n')
        self.f.write('write coor unit @fileu card \n\n')
        
    def minimization_mini(self):
        self.f.write('skip all excl bond \n')
        self.f.write('update inbfrq 0 1hbfrq 0 imgfrq 0 \n')
        self.f.write('mini sd nstep 100 \n')
        self.f.write('mini abnr nstep 1000000 tolg 0.00001 \n\n')
        
    def write_crd_mini(self, cutoff):
        self.f.write(f'open write unit 12 card name @nadir/cutoffdata/minim_after_fm_{cutoff:.2f}.crd\n')
        self.f.write('write coor card unit 12\n')
        self.f.write('close unit 12\n\n')
        
    def stream_str(self, cutoff):
        self.f.write('stream @topdir/na_enm_{0:.2f}.str \n\n'.format(cutoff))
        
    def stream_str_nma(self):
        self.f.write('stream @icstr \n\n')
        self.f.write('ic fill \n')
        self.f.write('open write unit @fileu card name @avgic \n')
        self.f.write('write ic   unit @fileu card resid \n')
        self.f.write('* Internal coordinate averages \n')
        self.f.write('* \n\n')
        self.f.write('close unit @fileu \n\n')
        
    def nma(self):
        self.f.write('calc nmode   ?natom * 3 \n\n')
        self.f.write('set nmodes   @nmode \n')
        self.f.write('set type     temp \n')
        self.f.write('set fluctu   20 \n')
        self.f.write('set vibu     22 \n\n')
        self.f.write('open write unit @fluctu card name @fluctic \n')
        self.f.write('open write unit @vibu   card name @vibic \n\n')
        self.f.write('vibran nmode @nmodes \n')
        self.f.write('    diag fini \n')
        self.f.write('    fluc ic @type @temp tfre 0.0 mode 7 thru @nmodes \n')
        self.f.write('    ic save \n')
        self.f.write('    ic write unit @fluctu resid \n')
        self.f.write('    * Internal coordinate fluctuation \n')
        self.f.write('    * \n')
        self.f.write('    write normal card mode 1 thru @nmodes unit @vibu \n')
        self.f.write('end \n\n')
        
    def read_traj(self):
        self.f.write('open unform read unit 21 name @dcddir/mode.@modeid.dcd \n')
        self.f.write('traj firstu 21 nunit 1 skip 1 \n\n')
     
    def icfluct(self):
        self.f.write('ic fill \n')
        self.f.write('ic dynam fluc firstu 21 nunit 1 begin 1 \n\n')
        
    def icavg(self, distance_average=False):
        self.f.write('ic fill \n')
        if distance_average:
            self.f.write('ic dynam avg firstu 21 nunit 1 begin 1 \n\n')
        
    def write_icfluct(self, cutoff):    
        self.f.write('open unit 30 write card name @rootfolder/@host/@typena/ic/mode.@modeid.{0:.2f}.ic \n'.format(cutoff))
        self.f.write('write ic card unit 30 \n')
        self.f.write('close unit 30 \n\n')
        
    def write_icavg(self, cutoff):    
        self.f.write('open unit 30 write card name @rootfolder/@host/@typena/ic/mode.@modeid.avg.{0:.2f}.ic \n'.format(cutoff))
        self.f.write('write ic card unit 30 \n')
        self.f.write('close unit 30 \n\n')
    
    def end(self):
        self.f.write('stop')
        self.f.close()
        
