import timeit
import re


# py fileproc.py out_files_corexome/out_sh/out_new_pos_mismatch_rs.txt coreexome_map
class NewPosMisM_06:

    def newposmims(self,brk=0):
        # variant position does not match db (b38) position(s)
        report = open('report_newposmims_' + self.ts + '.txt',"w")
        count = 0
        for line in self.inp.read():
            count+=1
            if brk and count == brk:
                break
            self.percent_comp(current=count,perc=10,total=brk)
            new_current = re.split('\tlookup=',line.rstrip())
            mid,brtab_hard,omtab_hard = re.split('\t',new_current[0])
            if not self.stillmain:
                print('omics db id %s not present (in consensus table) any longer' % (mid))
                continue
            brpos = self.checkbr_pos(mid)[1] #brtab_hard has single unmatching coord but some ids have multiple coords
            ompos = self.checkom_pos(mid)[1] #omtab_hard is unreliable, found to contain b37
            if len(new_current) == 1: # no look up available and probably not rsid
                nolkup = 0
                newtobr = [pos for pos in ompos if pos not in brpos and pos != '0:0']
                newtoom = [pos for pos in brpos if pos not in ompos and pos != '0:0']
                if newtoom:
                    report.write('No lookup available (re.split failed) for id %s. found position(s) new to omics: %s. Currently in omics: %s. Adding to omics positions table.\n' % (mid,','.join(newtoom),','.join(ompos)))
                    [self.addpos(mid,chrpos=cp,ds=self.tabname,build='38') for cp in newtoom]
                    nolkup += 1
                if newtobr:
                    print('No lookup available (re.split failed) for id %s. found position(s) new to br map table: %s. Currently in br map table: %s' % (mid,','.join(newtobr),','.join(brpos)))
                    nolkup += 1
                if not nolkup:
                    print('can remove this message after test stage: No lookup available for id %s. No novel coord found for either br table (%s) or omics db (%s).' % (mid,','.join(brpos),','.join(ompos)))
                continue
            b38 = self.getb38(new_current[1])
            if not b38:
                print('NO b38 FROM dbsnp AVAILABLE. %s. no action yet' % (mid))
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
                    report.write('br map table has coords for %s that do not match with dbsnp %s (%s), but also coords that do (%s). No action coded yet\n' % (mid,b38,','.join(brpos_mis),','.join(brpos_hit)))
                else:
                    report.write('br map table coords for %s (%s) do not match dbsnp coord (%s). Correcting the coordinates in map table\n' % (mid,','.join(brpos_mis),b38))
                    [self.mtab_change_pos(mid,oldpos=mis,newpos=b38) for mis in brpos_mis]
                actions += 1
            if ompos_mis:
                report.write('unchecked: omics b38 coord(s) for id %s (%s) do not match with the dbsnp coord (%s). Flagging entry in db\n' % (mid,','.join(ompos_mis),b38))
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
