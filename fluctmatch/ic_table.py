import numpy as np
import pandas as pd
from fluctmatch.fluctpair import FluctPair


class ICTable:
    def __init__(self, filename, skipheader=5, initial=False):
        self.filename = filename
        self.skipheader = skipheader
        if initial:
            self.ics, self.values = self.read_ic_initial()
        else:
            self.ics, self.values = self.read_ic()

    def read_ic(self):
        result = list()
        values = list()
        temp = np.genfromtxt(self.filename, skip_header=self.skipheader, dtype=str)
        for data in temp:
            result.append(FluctPair(data[3], data[7], float(data[17])))
            values.append(float(data[17]))
        return result, np.array(values)

    def read_ic_initial(self):
        result = list()
        values = list()
        temp = np.genfromtxt(self.filename, skip_header=self.skipheader, dtype=str)
        for data in temp:
            result.append(FluctPair(data[2], data[4], float(data[9])))
            values.append(float(data[9]))
        return result, np.array(values)

    def convert_to_df(self):
        d = {'I': list(), 'J': list(), 'r_IJ': list()}
        for pair, value in zip(self.ics, self.values):
            d['I'].append(pair.name1)
            d['J'].append(pair.name2)
            d['r_IJ'].append(value)
        df = pd.DataFrame(d)
        return df


class KBPair:
    def __init__(self, read_from_prm=False, filename=None, icavg=None, icfluct=None, rt=None):
        if read_from_prm:
            self.filename = filename
            self.d = self.read_prm()
        else:
            self.avg = icavg
            self.fluct = icfluct
            self.RT = rt
            self.d = self.read_ic()

    def read_prm(self):
        d = {'name1': list(), 'name2': list(), 'k': list(), 'b': list()}
        temp = np.genfromtxt(self.filename, skip_header=4, skip_footer=1, dtype=str)
        for data in temp:
            d['name1'].append(data[0])
            d['name2'].append(data[1])
            d['k'].append(float(data[2]))
            d['b'].append(float(data[3]))
        d['k'] = np.array(d['k'])
        d['b'] = np.array(d['b'])
        return d

    def read_ic(self):
        ks = self.RT / np.square(self.fluct.values)
        d = {'name1': list(), 'name2': list(), 'k': list(), 'b': list()}
        for pair, b, k in zip(self.avg.ics, self.avg.values, ks):
            d['name1'].append(pair.name1)
            d['name2'].append(pair.name2)
            d['b'].append(b)
            d['k'].append(k)
        d['k'] = np.array(d['k'])
        d['b'] = np.array(d['b'])
        return d

    def set_d_k(self, ks):
        self.d['k'] = ks

    def set_d_b(self, bs):
        self.d['b'] = bs

    def convert_to_df(self):
        return pd.DataFrame(self.d)



