import timeit
import re

# py fileproc.py out_files_corexome/out_sh/out_new_rs_rs.txt coreexome_map

class NewRs_07:
#map table rsid not known to db, map table non-rsid in db consensus instead
# related (newrsbyalt): map table rsid not known to db, map table non-rsid in db alt_ids instead
    def newrs(self,brk=0,start=0):
        report = open('report_newrs_' + self.ts + '_' + self.tabname + '.txt',"w")
        count = 0
        for line in self.inp.read():
            count+=1
            if brk and count == brk:
                break
            if start and count < start:
                continue
            self.percent_comp(current=count,perc=10,total=brk,start=start)
            new_current = re.split('\tlookup=',line.rstrip())
            mapid,linkid = re.split('\t',new_current[0])
            mapid = re.split('= ',mapid)[1]
            linkid = re.split('= ',linkid)[1]
            mgchk = self.mergecheck(new_current[1])
            merges = mgchk['merges']
            becomes = self.mergedinto(evenarr=merges,potbef=mapid)
            if becomes:
                print('map table id %s should be merged into %s. No action yet' % (mapid,'/'.join(becomes)))
                continue
#            if not self.stillmain(linkid):
#                print('link id %s is not in consensus table. Would have swapped it out with %s from map table' % (linkid,mapid))
#                continue
            b38 = self.getb38(new_current[1])
            getbrpos = self.checkbr_pos(linkid,b38)
            if getbrpos[0]:
                print('swapping main id %s out of omics consensus table and new id %s from map table in' % (linkid,mapid))
                self.swapout_main(swin=mapid,swout=linkid,ds=self.tabname)               
                continue
            print('no action ',mapid,linkid,b38,getbrpos)
