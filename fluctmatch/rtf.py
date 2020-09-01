import MDAnalysis

"""
https://www.charmm.org/charmm/documentation/by-version/c41b1/params/doc/io/#RTFFileformat
"""


class RTF:
    def __init__(self, host, type_na, crd, pairs, d_mass):
        self.host = host
        self.type_na = type_na
        self.crd = crd
        self.pairs = pairs
        self.u = MDAnalysis.Universe(self.crd, self.crd)
        self.atomnames = self.get_atomnames()
        self.title = self.make_title()
        self.mass = self.make_mass(d_mass)
        self.defa = ['DEFAULT FIRS NONE LAST NONE\n']
        self.res = ['RESI NA 0.00\n', 'GROUP\n']
        self.atoms = self.make_atoms()
        self.ics = self.make_ics()
        self.bonds = self.make_bonds()

    def get_atomnames(self):
        result = list()
        segid1 = self.u.select_atoms("segid STRAND1")
        for i, atom in enumerate(segid1):
            cgname = 'A{0}'.format(i+1)
            result.append(cgname)
        segid2 = self.u.select_atoms("segid STRAND2")
        for i, atom in enumerate(segid2):
            cgname = 'B{0}'.format(i+1)
            result.append(cgname)
        return result

    def make_title(self):
        return ['* {0}-{1}\n'.format(self.host, self.type_na), '*\n', '41  1\n']

    def make_mass(self, d_mass):
        result = list()
        idx = 1
        for atomname in self.atomnames:
            temp = 'MASS{0:6} {1:<6} {2:8.6f}\n'.format(idx, atomname, d_mass[atomname])
            result.append(temp)
            idx += 1
        return result

    def make_mass_all_1(self):
        result = list()
        idx = 1
        for atomname in self.atomnames:
            temp = 'MASS{0:6} {1:<6} {2:8.6f}\n'.format(idx, atomname, 1.000)
            result.append(temp)
            idx += 1
        self.mass = result

    def make_atoms(self):
        result = list()
        for atomname in self.atomnames:
            temp = 'ATOM {0:<5}{0:<8}0.00\n'.format(atomname)
            result.append(temp)
        return result

    def write_rtf(self, f_out):
        f = open(f_out, 'w')
        for text in self.title:
            f.write(text)
        f.write('\n')
        for text in self.mass:
            f.write(text)
        f.write('\n')
        for text in self.defa:
            f.write(text)
        for text in self.res:
            f.write(text)
        for text in self.atoms:
            f.write(text)
        for text in self.bonds:
            f.write(text)
        f.write('\n\n\nEND')
        #  for text in self.ics:
        #    f.write(text)
        f.close()

    def make_ics(self):
        result = list()
        n_pairs = len(self.pairs)
        last_idx = len(self.pairs) - 1
        if n_pairs % 2 == 1:
            odd = True
        else:
            odd = False
        idx_1 = 0
        idx_2 = 1
        for pair1 in self.pairs[::2]:
            if odd and idx_1 == last_idx:
                pair2 = self.pairs[0]
            else:
                pair2 = self.pairs[idx_2]
            temp1 = '{0} {1} {2} {3}'.format(pair1.name1, pair1.name2, pair2.name1,
                                             pair2.name2)
            temp2 = 'IC {0} 1.0000 0.00 180.00 0.00 0.0000\n'.format(temp1)
            result.append(temp2)
            idx_1 += 2
            idx_2 += 2
        return result

    def make_bonds(self):
        result = list()
        for pair in self.pairs:
            temp = 'BOND {0:5} {1:<5}\n'.format(pair.name1, pair.name2)
            result.append(temp)
        return result

