from .types import Types
from datetime import datetime
import sys
import os
import re
from utils.queryfile import QueryFile, NormFile

class ProcFile(Types):

    def __init__(self,fname,tabname,omics,br):
        self.inp = NormFile(fname)
        self.ts = datetime.now().strftime("%d%b_%Hhr")
        self.omics = omics
        self.br = br
        self.tabname = tabname
        print("line count",self.inp.row_count)
        if 'map_new_alt_rs' in self.inp.bfile:
            if 'out' in os.path.basename(self.inp.bfile):
                self.newaltrs()
            else:
                sys.exit("for file \"%s\" there should be an equivalent with dbsnp look up (%s %s -> %s)" % (self.inp.bfile,'map_new_alt_rs.sh',self.inp.bfile,'out_map_new_alt_rs.txt'))

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
            return {'withdrawn':False,'merges':[]}
        if 'withdrawn' in rside:
            return {'withdrawn':True}
        mg = re.findall("(rs[0-9]+)",rside)
        if (len(mg) % 2) != 0:
            raise Exception("mergecheck(col): expecting paired (even) rsids in merge list: ", col)
        return {'withdrawn':False,'merges':mg}

    def frominto(self,evenarr,bfor,aftr):
        #print(evenarr)
        nxt = []
        correct = False
        after_match = False
        for ind,rs in enumerate(evenarr):
            ind += 1
            if (ind % 2) != 0:
                if rs == bfor:
                    correct = True
                    nxt.append(evenarr[ind])
        if correct:
            if aftr in nxt:
                after_match = True
        return {'correct':correct,'after_match':after_match,'alt_after':nxt}

    def getb38(self,container):
        rslt = re.search('(?:b38=)(\S+)',container)
        return rslt.group(1)

    def add_alt(self,alt,main):
        pass

    def mtab_change_id(self,xsting,chngto):
        pass

    def swapout_main(self,swin,swout):
        pass

    def checkbr_pos(self,mid,chrpos):
        if 'ukbbaffy_v2_1_map' in self.tabname:
            rscol = 'chipid'
        else:
            rscol = 'dbsnpid'
        q = 'SELECT chr FROM ' + self.tabname + ' WHERE snp = %s'
        val = (mid,)
        self.br.execute(q,val)
        res = self.br.fetchall()
        if self.br.rowcount == 0 and 'rs' in mid: # or len(res)
            q = 'SELECT chr FROM ' + self.tabname + ' WHERE ' + rscol + ' REGEXP %s'
            val = (mid+'[[:>:]]',)
            self.br.execute(q,val)
            res = self.br.fetchall()
        res = [s[0] for s in res]
        if chrpos in res:
            return True
        return False

    def checkom_pos(self,mid,chrpos):
        q = 'SELECT chr,pos FROM positions WHERE build = %s AND id = %s'
        val = ('38',mid)
        self.omics.execute(q,val)
        res = self.omics.fetchall()
        res = [s[0] + ':'+ str(s[1]) for s in res]
        if chrpos in res:
            return True
        return False

    def extra_map(self,newid,linkid,chrpos,datasource,chosen):
        pass
