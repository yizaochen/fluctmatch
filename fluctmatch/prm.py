"""
https://www.charmm.org/charmm/documentation/by-version/c41b1/params/doc/io/#RTFFileformat
"""


class PRM:
    def __init__(self, host, type_na, kbpairs, iternum=0):
        self.host = host
        self.type_na = type_na
        self.name1s = kbpairs.d['name1']
        self.name2s = kbpairs.d['name2']
        self.bs = kbpairs.d['b']
        self.ks = kbpairs.d['k']
        self.n_iter = iternum
        self.title = self.make_title()
        self.bonds = self.make_bonds()

    def make_title(self):
        return ['* {0}-{1} iternum={2}\n'.format(self.host, self.type_na, self.n_iter),
                '*\n\n']

    def write_prm(self, f_out):
        f = open(f_out, 'w')
        for text in self.title:
            f.write(text)
        for text in self.bonds:
            f.write(text)
        f.write('\nEND')
        f.close()

    def make_bonds(self):
        result = ['BONDS\n']
        for name1, name2, k, b in zip(self.name1s, self.name2s, self.ks, self.bs):
            temp = '{0:<7}{1:<6}{2:9.4f}  {3:8.4f}\n'.format(name1, name2, k, b)
            result.append(temp)
        return result

