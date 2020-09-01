import MDAnalysis


class ICSTR:
    def __init__(self, host, type_na, crd, pairs, atomid_map):
        self.host = host
        self.type_na = type_na
        self.crd = crd
        self.pairs = pairs
        self.atomid_map = atomid_map
        self.u = MDAnalysis.Universe(self.crd, self.crd)
        self.title = self.make_title()
        self.ics = self.make_ics()

    def make_title(self):
        return ['* {0}-{1}\n'.format(self.host, self.type_na)]

    def make_ics(self):
        i = 1
        result = ['IC EDIT\n']
        for pair in self.pairs:
            if i % 300 == 0:
                result.append('END\n')
                result.append('IC EDIT\n')
            id_1 = self.atomid_map[pair.name1]
            id_2 = self.atomid_map[pair.name2]
            temp = 'DIST BYNUM{0:10d} BYNUM{1:10d}\n'.format(id_1, id_2)
            result.append(temp)
            i += 1
        return result

    def write_ic_str(self, f_out):
        f = open(f_out, 'w')
        for text in self.title:
            f.write(text)
        for text in self.ics:
            f.write(text)
        f.write('END')
        f.close()


