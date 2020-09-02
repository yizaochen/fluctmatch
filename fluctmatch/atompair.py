import numpy as np

class AtomPair:
    def __init__(self, name1, name2, selection1, selection2, k=10., b=5.):
        self.name1 = name1
        self.name2 = name2
        self.selection1 = selection1
        self.selection2 = selection2
        self.atom1 = None  # MDAnalysis Atom Obejct
        self.atom2 = None  # MDAnalysis Atom Obejct
        self.distances = list()
        self.k = k
        self.b = b
        self.pairname = "{0}-{1}".format(self.name1, self.name2)

    def set_atom1(self, atom):
        self.atom1 = atom

    def set_atom2(self, atom):
        self.atom2 = atom

    def get_distance(self):
        return np.linalg.norm(self.atom1.positions[0] - self.atom2.positions[0])

    def append_distances(self):
        self.distances.append(np.linalg.norm(self.atom1.positions[0] - self.atom2.positions[0]))

    def __repr__(self):
        return "{0}-{1}".format(self.name1, self.name2)

    def __hash__(self):
        return hash(self.name1) + hash(self.name2)

    def __eq__(self, other):
        if self.name1 == other.name1 and self.name2 == other.name2:
            return True
        elif self.name1 == other.name2 and self.name2 == other.name1:
            return True
        else:
            return False