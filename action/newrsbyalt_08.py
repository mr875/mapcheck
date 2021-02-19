import timeit
import re

# py fileproc.py out_files_corexome/out_sh/out_new_rs_byalt_rs.txt coreexome_map
# py fileproc.py out_files_humanexome/out_sh/out_new_rs_byalt_rs.txt humanexome_map
# py fileproc.py out_files_infimun/out_sh/out_new_rs_byalt_rs.txt infiniumimmunoarray_map
# py fileproc.py out_files_msex/out_sh/out_new_rs_byalt_rs.txt msexome_map

class NewRsbyAlt_08:
# map table rsid not known to db, map table non-rsid in db alt_ids instead
    def newrsbyalt(self,brk=0,start=0):
        report = open('report_newrsbyalt_' + self.ts + '_' + self.tabname + '.txt',"w")
        count = 0
        for line in self.inp.read():
            count+=1
            if brk and count == brk:
                break
            if start and count < start:
                continue
            self.percent_comp(current=count,perc=10,total=brk,start=start)
            new_current = re.split('\tlookup= ',line.rstrip())
            mapid,linkid,midom = re.split('\t',new_current[0])
            mapid = re.split('= ',mapid)[1]
            linkid = re.split('= ',linkid)[1]
            midom = re.split('= ',midom)[1]
            mgchk = self.mergecheck(new_current[1])
            merges = mgchk['merges']
            #print('map id = %s\tlink id = %s\tmain omics id = %s' % (mapid,linkid,midom))
            becomes = self.mergedinto(evenarr=merges,potbef=mapid)
            stillmain = self.stillmain(midom)
            mrgid = None
            if becomes:
                if len(becomes) > 1: # probably not necessary
                    print('map table id %s is merged into multiple (%s) according to dbsnp lookup. No action yet' % (mapid,'/'.join(becomes)))
                    continue
                mrgid = becomes[0]
                report.write('map table id %s merged to %s and linked to omics %s via %s. correcting map table (%s -> %s)\n' % (mapid,mrgid,midom,linkid,mapid,mrgid))
                self.mtab_change_id(xsting=mapid,chngto=mrgid)
                mrgid_knwn = self.checkomid(cid=mrgid) # is the new mrgid unknown as the mapid was?
                if not mrgid_knwn[0] and not mrgid_knwn[1]:
                    print('new id %s is also unknown to omics so adding old map id %s as alternative id to %s (if still present) and continuing\n' % (mrgid,mapid,midom))
                    if stillmain:
                        self.add_alt(alt=mapid,main=midom,ds=self.tabname)
                else:
                    report.write('new merged id %s is already known as consensus id or alternative id or both (%s). No action\n' % (mrgid,'-'.join([str(i) for i in mrgid_knwn])))
                    continue
            if not stillmain:
                report.write('main id %s (linked by %s) is not in consensus table anymore. Would have swapped it out with %s from map table (or merged %s)\n' % (midom,linkid,mapid,mrgid))
                continue
            b38 = self.getb38(new_current[1])
            if not b38:
                if not mrgid:
                    report.write('b38 position of map table id %s (omics %s, linked by %s) not available (may be withdrawn). Will add to omics as alternative id instead\n' % (mapid,midom,linkid))
                    self.add_alt(alt=mapid,main=midom,ds=self.tabname)
                else:
                    print('b38 position of new merge id %s (omics %s, linked by %s) not available (may be withdrawn). Will add to omics as alternative id instead\n' % (mrgid,midom,linkid))
                    self.add_alt(alt=mrgid,main=midom,ds='dbsnp')
                continue
            contig = False
            if '_' in b38:
                report.write('map table id %s (omics %s, linked by %s) dbsnp look up has a contig id instead of regular chromosome identifier: %s. Will swap %s (or merge %s) into omics for %s but will not correct/update positions\n' % (mapid,midom,linkid,b38,mapid,mrgid,linkid))
                contig = True
            getbrpos = self.checkbr_pos(mapid,b38)
            getompos = self.checkom_pos(midom,chrpos=b38,build='38')
            if not contig:
                if not getbrpos[0]:
                    oldbrpos = [op for op in getbrpos[1] if op != '0:0']
                    for op in oldbrpos:
                        report.write('map table id %s (or merged %s) to be swapped in to consensus table in place of %s (linked by %s). %s in map table does not match dbsnp coordinate (%s vs %s). Correcting position in map table\n' % (mapid,mrgid,midom,linkid,mapid,op,b38))
                        self.mtab_change_pos(anid=mapid,oldpos=op,newpos=b38)
                if not getompos[0]:
                    report.write('map table id %s (or merged %s) to be swapped in to consensus table in place of %s (linked by %s). %s in omics does not match dbsnp coordinate (%s vs %s). Adding dbsnp coordinate %s to omics\n' % (mapid,mrgid,midom,linkid,midom,b38,','.join(getompos[1]),b38))
                    badom = [bp for bp in getompos[1] if bp != '0:0']
                    self.addpos(midom,chrpos=b38,ds='dbsnp',build='38')
                    for bp in badom:
                        print('flagging position %s in omics db because we have a dbsnp sourced coordinate %s for this variant (%s/%s in omics, %s/%s in map table)\n' % (bp,b38,midom,linkid,mapid,mrgid))
                        self.pos_flag(midom,chrpos=bp,fl=-5)
            if not mrgid:
                report.write('swapping main id %s (linked by %s) out of omics consensus table and new id %s from map table in\n' % (midom,linkid,mapid))
                self.swapout_main(swin=mapid,swout=midom,ds=self.tabname) 
            else:
                print('swapping main id %s (linked by %s) out of omics consensus table and new merged id %s from map table (originally %s) in\n' % (midom,linkid,mrgid,mapid))
                self.swapout_main(swin=mrgid,swout=midom,ds='dbsnp') 
        print('Done. %s seconds' % (int(timeit.default_timer() - self.grandstart)))
