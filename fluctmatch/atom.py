class Atom:
    def __init__(self, data, segid_exist):
        self.record_name = data[0]
        self.serial = int(data[1])
        self.name = data[2]
        self.resname = data[3]
        self.resid = int(data[4])
        self.x = float(data[5])
        self.y = float(data[6])
        self.z = float(data[7])
        self.occupancy = float(data[8])
        self.tempFactor = float(data[9])
        if segid_exist:
            self.segid = data[10]

    def set_resid(self, resid):
        self.resid = int(resid)

    def set_resname(self, resname):
        self.resname = resname

    def set_atomid(self, atomid):
        self.serial = int(atomid)

    def set_tempFactor(self, tempFactor):
        self.tempFactor = tempFactor
        
    def set_segid(self, segid):
        self.segid = segid

    def get_format_str_pdb(self):
        if len(self.name) < 4:
            str1 = f'{self.record_name:4s}{self.serial:>7d}  {self.name:<3s}'
        else:
            str1 = f'{self.record_name:4s}{self.serial:>7d} {self.name:4s}'
        str2 = f' {self.resname:<4s}  {self.resid:>3d}'
        str3 = f'{self.x:>12.3f} {self.y:>7.3f} {self.z:>7.3f}'
        str4 = f'  1.00 {self.tempFactor:>6.2f}      '
        return str1 + str2 + str3 + str4
