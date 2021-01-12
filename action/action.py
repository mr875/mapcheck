#from utils.connect import DBConnect
from utils.queryfile import QueryFile, NormFile

class ProcFile:

    def __init__(self,fname):
        self.inp = NormFile(fname)
        print("line count",self.inp.row_count)
        if 'map_new_alt_rs' in self.inp.bfile:
            self.newaltrs()

    def newaltrs(self):
        count = 0
        for line in self.inp.readls():
            count+=1
            print(line[0])
            if count == 10:
                break

