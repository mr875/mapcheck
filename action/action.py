#from utils.connect import DBConnect
import sys
import os
import re
from utils.queryfile import QueryFile, NormFile

class ProcFile:

    def __init__(self,fname):
        self.inp = NormFile(fname)
        print("line count",self.inp.row_count)
        if 'map_new_alt_rs' in self.inp.bfile:
            if 'out' in os.path.basename(self.inp.bfile):
                self.newaltrs()
            else:
                sys.exit("for file \"%s\" there should be an equivalent with dbsnp look up (%s %s -> %s)" % (self.inp.bfile,'map_new_alt_rs.sh',self.inp.bfile,'out_map_new_alt_rs.txt'))

    def newaltrs(self):
        count = 0
        for line in self.inp.read():
            count+=1
            print(line)
            new_current = re.split('Current',line)
            newrs = self.grabrsinp(new_current[0])
            currs = self.grabrsinp(new_current[1])
            newmchk = self.mergecheck(new_current[0])
            curmchk = self.mergecheck(new_current[1])
            #self.frominto(newmchk,newrs)
            if count == 6:
                break

    def grabrsinp(self,col):
        rs = re.search("(?:\()(rs[0-9]+)(?:\))", col)
        if rs:
            return rs.group(1)
        else:
            raise Exception("grabrsinp(col): can't find rs id in parenthesis: %s" % (col))

    def mergecheck(self,col):
        rside = re.search("(?:[mM]erges?:)(.+)",col)
        rside = rside.group(1)
        if 'rs' not in rside:
            return None
        if 'withdrawn' in rside:
            return {'withdrawn':True}
        mg = re.findall("(rs[0-9]+)",rside)
        if (len(mg) % 2) != 0:
            raise Exception("mergecheck(col): expecting paired (even) rsids in merge list: ", col)
        return {'withdrawn':False,'merges':mg}

    def frominto(self,evenarr,qrs):
        print(evenarr)
        for ind,rs in enumerate(evenarr):
            ind += 1
