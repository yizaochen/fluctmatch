import numpy as np
from fluctmatch import atom

class PDBReader:
    """
    Reference: https://www.mdanalysis.org/docs/_modules/MDAnalysis/coordinates/PDB.html#PDBReader

    Read a PDB file and output a list contains a lot of atom objects
    """
    def __init__(self, fname, skip_header=1, skip_footer=1, segid_exist=False):
        self.fname = fname
        self.skip_header = skip_header
        self.skip_footer = skip_footer
        self.segid_exist = segid_exist
        self.atomgroup = self._read_pdb()

    def get_atomgroup(self):
        return self.atomgroup

    def _read_pdb(self):
        lines = np.genfromtxt(self.fname, skip_header=self.skip_header, skip_footer=self.skip_footer, dtype=str)
        atomgroup = list()
        for line in lines:
            atomgroup.append(atom.Atom(line, self.segid_exist))
        return atomgroup


class PDBWriter:
    def __init__(self, fname, atomgroup):
        self.fname = fname
        self.atomgroup = atomgroup

    def write_pdb(self):
        f = open(self.fname, 'w')
        f.write('REMARK\n')
        for atomobj in self.atomgroup:
            strout = atomobj.get_format_str_pdb()
            f.write(f'{strout}\n')
        f.write('END')
        f.close()
        print(f'Write PDB: {self.fname}')

