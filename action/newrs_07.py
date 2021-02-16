import timeit
import re

# py fileproc.py out_files_corexome/out_sh/out_new_rs_rs.txt coreexome_map

class NewRs_07:
# map table rsid not known to db, map table non-rsid in db consensus instead
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
            if not self.stillmain(linkid):
                print('link id %s is not in consensus table. Would have swapped it out with %s from map table' % (linkid,mapid))
                continue
            b38 = self.getb38(new_current[1])
            if not b38:
                print('b38 position of map table id %s (omics %s) not available (may be withdrawn)' % (mapid,linkid))
                continue
            contig = False
            if '_' in b38:
                print('map table id %s (omics %s) dbsnp look up has a contig id instead of regular chromosome identifier: %s. Will swap %s into omics for %s but will not correct/update positions' % (mapid,linkid,b38,mapid,linkid))
                contig = True
            getbrpos = self.checkbr_pos(mapid,b38)
            getompos = self.checkom_pos(linkid,chrpos=b38,build='38')
            if not contig:
                if not getbrpos[0]:
                    oldbrpos = [op for op in getbrpos[1] if op != '0:0']
                    for op in oldbrpos:
                        report.write('map table id %s to be swapped in to consensus table in place of %s. %s in map table does not match dbsnp coordinate (%s vs %s). Correcting position in map table\n' % (mapid,linkid,mapid,op,b38))
                        self.mtab_change_pos(anid=mapid,oldpos=op,newpos=b38)
                if not getompos[0]:
                    report.write('map table id %s to be swapped in to consensus table in place of %s. %s in omics does not match dbsnp coordinate (%s vs %s). Adding dbsnp coordinate %s to omics\n' % (mapid,linkid,linkid,b38,','.join(getompos[1]),b38))
                    badom = [bp for bp in getompos[1] if bp != '0:0']
                    self.addpos(linkid,chrpos=b38,ds='dbsnp',build='38')
                    for bp in badom:
                        print('flagging position %s in omics db because we have a dbsnp sourced coordinate %s for this variant (%s/%s)' % (bp,b38,linkid,mapid))
                        self.pos_flag(linkid,chrpos=bp,fl=-5)
            report.write('swapping main id %s out of omics consensus table and new id %s from map table in\n' % (linkid,mapid))
            self.swapout_main(swin=mapid,swout=linkid,ds=self.tabname)               
        print('Done. %s seconds' % (int(timeit.default_timer() - self.grandstart)))
