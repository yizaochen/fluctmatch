from fluctmatch.atompair import AtomPair

class FluctPair(AtomPair):
    def __init__(self, name1, name2, fluct_value):
        self.name1 = name1
        self.name2 = name2
        self.value = fluct_value