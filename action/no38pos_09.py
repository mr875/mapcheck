import timeit
import re

# py fileproc.py out_files_corexome/out_sh/out_no38pos_rs.txt coreexome_map
# py fileproc.py out_files_msex/out_sh/out_no38pos_rs.txt msexome_map
# py fileproc.py out_files_omni/out_sh/out_no38pos_rs.txt omniexpress_map
# py fileproc.py out_files_humanexome/out_sh/out_no38pos_rs.txt humanexome_map

# py fileproc.py out_files_infimun/out_sh/out_no38pos_rs.txt infiniumimmunoarray_map
# py fileproc.py out_files_omni21/out_sh/out_no38pos_rs.txt omniexpress_v2_1_map
# py fileproc.py out_files_ukbb2_1/out_sh/out_no38pos_rs.txt ukbbaffy_v2_1_map

class No38Pos_09:
# variant found in db but no b38 position available
    def no38pos(self,brk=0,start=0):
        report = open('report_no38pos_' + self.ts + '_' + self.tabname + '.txt',"w")
        count = 0
        for line in self.inp.read():
            count+=1
            if brk and count == brk:
                break
            if start and count < start:
                continue
            self.percent_comp(current=count,perc=10,total=brk,start=start)
            new_current = re.split('\tlookup:',line.rstrip())
            mid,brtab_pos,brorig = re.split('\t',new_current[0])
            brorig = re.split('=',brorig)[1]
            if not self.stillmain(mid):
                report.write('omics db id %s not present (in consensus table) any longer. skipping\n' % (mid))
                continue
            ompos = self.checkom_pos(mid=mid,build='38')[1] 
            b38 = None
            becomes = []
            getbrpos = self.checkbr_pos(mid)
            brid = mid
            if not getbrpos and brorig != mid:
                getbrpos =self.checkbr_pos(brorig)
                brid = brorig
                report.write('mid %s maps to br map table id %s\n' % (mid,brid))
            else:
                brorig = None
            if len(new_current) > 1: # is dbsnp look up available 
                b38 = self.getb38(new_current[1])
                merges = self.mergecheck(new_current[1])['merges']
                if merges:
                    becomes = self.mergedinto(evenarr=merges,potbef=mid)
                    # ^^ check later. will contain new merge id if anything
            if not b38: # look up not available (and probably no rs id) OR b38 not available
                ompos37 = self.checkom_pos(mid=mid,build='37')[1] 
                if brtab_pos in ompos37:
                    report.write('mid: %s. No dbsnp lookup available. map table pos (from file) %s is present in omics as build 37 however (%s). No action because map table may not have b38 coordinates\n' % (mid,brtab_pos,' & '.join(ompos37)))
                else:
                    if brtab_pos in ompos:
                        report.write('mid %s. No dbsnp lookup available. map table pos (from file) %s was due to be added to omics under build 38 but it appears to exist already (%s). No action\n' % (mid,brtab_pos,' & '.join(ompos)))
                    else:
                        if brtab_pos != '0:0':
                            report.write('mid: %s. No dbsnp lookup available. map table pos (from file) %s to be added to omics under build 38\n' % (mid,brtab_pos))
                            self.addpos(mid=mid,chrpos=brtab_pos,ds=self.tabname,build='38')
                        else:
                            report.write('mid: %s. No dbsnp lookup available. map table pos (from file) is %s, and because it is 0:0 it will not be added to omics\n' % (mid,brtab_pos))
            else:
                #print('mid: %s. Lookup available: dbsnp b38 is %s' % (mid,b38))
                brpos = getbrpos[1] #brtab_hard has single unmatching coord but some ids have multiple coords
                ompos_mis = [pos for pos in ompos if pos != b38 and pos != '0:0']
                ompos_hit = [pos for pos in ompos if pos == b38]
                brpos_mis = [pos for pos in brpos if pos != b38 and pos != '0:0']
                brpos_hit = [pos for pos in brpos if pos == b38]
                posaction = 0
                if brpos_mis:
                    if brpos_hit:
                        report.write('br map table has coords for %s that do not match with dbsnp %s (%s), but also coords that do (%s). Correcting the non matching coordinates in map table\n' % (brid,b38,','.join(brpos_mis),','.join(brpos_hit)))
                    else:
                        report.write('br map table coords for %s (%s) do not match dbsnp coord (%s). Correcting the coordinates in map table\n' % (brid,','.join(brpos_mis),b38))
                    for mis in brpos_mis:
                        self.mtab_change_pos(brid,oldpos=mis,newpos=b38)
                    posaction += 1
                if ompos_mis:
                    print('omics b38 coord(s) for id %s (%s) do not match with the dbsnp coord (%s). Flagging entry in db\n' % (mid,','.join(ompos_mis),b38))
                    for mis in ompos_mis:
                        self.pos_flag(mid,mis,fl=-5)
                    posaction += 1
                if not ompos_hit:
                    report.write('dbsnp coord for id %s (%s) is not in omics db. Adding it.\n' % (mid,b38))
                    self.addpos(mid,chrpos=b38,ds='dbsnp',build='38')
                    posaction += 1
                if posaction == 0:
                    report.write('no position-based action done on %s' % (mid))
            if becomes:
                self.dbsnp_merge(becomes,mid,brid,report)
        print('Done. %s seconds' % (int(timeit.default_timer() - self.grandstart)))
        report.close()


    def dbsnp_merge(self,becomes,mid,brid,report):
        if len(becomes) > 1: # probably not necessary
            print('map table id %s is merged into multiple (%s) according to dbsnp lookup. No action yet' % (brid,'/'.join(becomes)))
            return
        mrgid = becomes[0]
        report.write('map table id %s merged to %s. correcting map table (%s -> %s)\n' % (brid,mrgid,brid,mrgid))
        self.mtab_change_id(xsting=brid,chngto=mrgid)
        mrgid_knwn = self.checkomid(cid=mrgid) 
        if mrgid_knwn[0]:
            report.write('mid %s is merged in dbsnp to %s, but %s already exists in consensus table. Adding as alternative id instead\n' % (mid,mrgid,mrgid))
            self.add_alt(alt=mrgid,main=mid,ds='dbsnp')
        else:
            print('mid %s is merged in dbsnp to %s. Swapping out %s and swapping in %s\n' % (mid,mrgid,mid,mrgid))
            self.swapout_main(swin=mrgid,swout=mid,ds='dbsnp')               
