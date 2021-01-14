#from utils.connect import DBConnect
from datetime import datetime
import sys
import os
import re
from utils.queryfile import QueryFile, NormFile

class ProcFile:

    def __init__(self,fname):
        self.inp = NormFile(fname)
        self.ts = datetime.now().strftime("%d%b_%Hhr")
        print("line count",self.inp.row_count)
        if 'map_new_alt_rs' in self.inp.bfile:
            if 'out' in os.path.basename(self.inp.bfile):
                self.newaltrs()
            else:
                sys.exit("for file \"%s\" there should be an equivalent with dbsnp look up (%s %s -> %s)" % (self.inp.bfile,'map_new_alt_rs.sh',self.inp.bfile,'out_map_new_alt_rs.txt'))

    def newaltrs(self):
        #rs id not found in db but db already has another rs id for that variant
        report = open('report_newaltr_' + self.ts + '.txt',"w")
        count = 0
        for line in self.inp.read():
            count+=1
            if count == 7:
                break
            #print(line)
            new_current = re.split('Current',line)
            newrs = self.grabrsinp(new_current[0])
            currs = self.grabrsinp(new_current[1])
            newmchk = self.mergecheck(new_current[0])
            curmchk = self.mergecheck(new_current[1])
            if newmchk['withdrawn']:
                if not curmchk['withdrawn']:
                    report.write('map table %s is withdrawn from dbsnp, use omics db %s instead\n' % (newrs,currs))
                    self.add_alt(alt=newrs,main=currs)
                    self.mtab_change_id(xsting=newrs,chngto=currs)
                else:
                    report.write('map table %s AND omics db %s are withdrawn from dbsnp\n' % (newrs,currs))
                continue   
            if curmchk['withdrawn']:
                report.write('omics db %s is withdrawn. map table %s to be used in omics db\n' % (currs,newrs))
                self.swapout_main(swin=newrs,swout=currs)
                continue
            mergelist = newmchk['merges'] + curmchk['merges']
            newtocurr = self.frominto(evenarr=mergelist,bfor=newrs,aftr=currs)
            if newtocurr['correct']:
                report.write("map table %s merged into omics db %s, use %s\n" % (newrs,currs,currs))
                self.add_alt(alt=newrs,main=currs)
                self.mtab_change_id(xsting=newrs,chngto=currs)
                continue
            currtonew = self.frominto(evenarr=mergelist,bfor=currs,aftr=newrs)
            if currtonew['correct']:
                report.write("omics db %s merged into map table %s, omics db to take %s\n" % (currs,newrs,newrs))
                self.swapout_main(swin=newrs,swout=currs)
                continue
            newrsb38 = self.getb38(new_current[0]) 
            print("no action coded for map table %s and omics db %s" % (newrs,currs))
        report.close()

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
        print(rslt.group(1))

    def add_alt(self,alt,main):
        pass

    def mtab_change_id(self,xsting,chngto):
        pass

    def swapout_main(self,swin,swout):
        pass
