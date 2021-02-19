import timeit
import re


# py fileproc.py out_files_corexome/out_sh/out_new_pos_mismatch_rs.txt coreexome_map
# py fileproc.py out_files_humanexome/out_sh/out_new_pos_mismatch_rs.txt humanexome_map
# py fileproc.py out_files_msex/out_sh/out_new_pos_mismatch_rs.txt msexome_map
# py fileproc.py out_files_omni/out_sh/out_new_pos_mismatch_rs.txt omniexpress_map
# py fileproc.py out_files_omni21/out_sh/out_new_pos_mismatch_rs.txt omniexpress_v2_1_map
# py fileproc.py out_files_ukbb2_1/out_sh/out_new_pos_mismatch_rs.txt ukbbaffy_v2_1_map 
# py fileproc.py out_files_infimun/out_sh/out_new_pos_mismatch_rs.txt infiniumimmunoarray_map 

class NewPosMisM_06:

    def newposmims(self,brk=0,start=0):
        # variant position does not match db (b38) position(s)
        report = open('report_newposmims_' + self.ts + '_' + self.tabname + '.txt',"w")
        count = 0
        for line in self.inp.read():
            count+=1
            if brk and count == brk:
                break
            if start and count < start:
                continue
            self.percent_comp(current=count,perc=10,total=brk,start=start)
            new_current = re.split('\tlookup=',line.rstrip())
            mid,brtab_hard,omtab_hard,brorig = re.split('\t',new_current[0])
            brorig = brorig.split('=')[1]
            if not self.stillmain:
                print('omics db id %s not present (in consensus table) any longer' % (mid))
                continue
            getbrpos = self.checkbr_pos(mid)
            brid = mid
            if not getbrpos and brorig != mid:
                getbrpos =self.checkbr_pos(brorig)
                brid = brorig
                report.write('mid %s maps to br map table id %s\n' % (mid,brid))
            else:
                brorig = None
            brpos = getbrpos[1] #brtab_hard has single unmatching coord but some ids have multiple coords
            ompos = self.checkomics_pos(mid)[0] #omtab_hard is unreliable, found to contain b37
            ompos = list(dict.fromkeys(ompos))
            b38 = None
            explained = 'No lookup available (re.split failed) '
            if len(new_current) > 1: # look up available 
                b38 = self.getb38(new_current[1])
                if not b38:
                    explained = 'No b38 from dbsnp (possibly withdrawn) '
            if not b38: # look up not available (and probably no rs id) OR b38 not available
                nolkup = 0
                newtobr = [pos for pos in ompos if pos not in brpos and pos != '0:0']
                newtoom = [pos for pos in brpos if pos not in ompos and pos != '0:0']
                if newtoom:
                    for cp in newtoom:
                        if self.checkom_pos(mid,chrpos=cp,build='37')[0]:
                            report.write(explained + 'for id %s. found position(s) new to omics: %s. Currently in omics (b38): %s. BUT position %s is in omics under b37. the position will not be added to omics again because it can not be validated by dbsnp look up.\n' % (mid,cp,','.join(ompos),cp))
                        else:
                            report.write(explained + 'for id %s. found position(s) new to omics: %s. Currently in omics: %s. Adding to omics positions table.\n' % (mid,cp,','.join(ompos)))
                            self.addpos(mid,chrpos=cp,ds=self.tabname,build='38')
                    nolkup += 1
                if newtobr:
                    newtobr_unflagged = [ntbr for ntbr in newtobr if not self.checkom_flag(mid,chrpos=ntbr)[0]]
                    oldbrpos = [obrp for obrp in brpos if obrp != '0:0']
                    if len(newtobr_unflagged) == 1:
                        ntbr = newtobr_unflagged[0]
                        report.write(explained + 'for id %s. found position new to br map table: %s. Currently in br map table: %s. No flag in omics so switching the position into br map table\n' % (mid,ntbr,','.join(brpos)))
                        for oldp in oldbrpos:
                            self.mtab_change_pos(brid,oldpos=oldp,newpos=ntbr)
                    else:
                        print(explained + 'for id %s. found position new to br map table: %s. Currently in br map table: %s. No change made in br map table due to omics entry being flagged or ambiguous (more than 1) new positions\n' % (mid,','.join(newtobr),','.join(brpos)))
                    nolkup += 1
                if not nolkup:
                    print('can remove this message after test stage: No lookup available for id %s. No novel coord found for either br table (%s) or omics db (%s).' % (mid,','.join(brpos),','.join(ompos)))
                continue
            mgchk = self.mergecheck(new_current[1])
            merges = mgchk['merges']
            brpos_mis = [pos for pos in brpos if pos != b38 and pos != '0:0']
            brpos_hit = [pos for pos in brpos if pos == b38]
            ompos_mis = [pos for pos in ompos if pos != b38 and pos != '0:0']
            ompos_hit = [pos for pos in ompos if pos == b38]
            actions = 0
            if brpos_mis:
                if brpos_hit:
                    report.write('br map table has coords for %s that do not match with dbsnp %s (%s), but also coords that do (%s). Correcting the non matching coordinates in map table\n' % (brid,b38,','.join(brpos_mis),','.join(brpos_hit)))
                else:
                    report.write('br map table coords for %s (%s) do not match dbsnp coord (%s). Correcting the coordinates in map table\n' % (brid,','.join(brpos_mis),b38))
                for mis in brpos_mis:
                    self.mtab_change_pos(brid,oldpos=mis,newpos=b38)
                actions += 1
            if ompos_mis:
                report.write('omics b38 coord(s) for id %s (%s) do not match with the dbsnp coord (%s). Flagging entry in db\n' % (mid,','.join(ompos_mis),b38))
                [self.pos_flag(mid,chrpos,fl=-5) for chrpos in ompos_mis]
                actions += 1
            if not ompos_hit:
                report.write('dbsnp coord for id %s (%s) is not in omics db. Adding it.\n' % (mid,b38))
                self.addpos(mid,chrpos=b38,ds='dbsnp',build='38')
                actions += 1
            if actions == 0:
                print('no position-based action done on %s' % (mid))
            becomes = self.mergedinto(evenarr=merges,potbef=mid)
            if becomes:
                print('existing main id %s should be merged into %s. Not action yet' % (mid,'/'.join(becomes)))
        print('Done. %s seconds' % (int(timeit.default_timer() - self.grandstart)))
        report.close()
