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
            if len(new_current) == 1: # no look up available and probably not rsid
                continue
            b38 = self.getb38(new_current[1])
            mgchk = self.mergecheck(new_current[1])
            merges = mgchk['merges']
            brpos = self.checkbr_pos(mid)[1] #brtab_hard has single unmatching coord but some ids have multiple coords
            ompos = self.checkom_pos(mid)[1] #omtab_hard is unreliable, found to contain b37
            brpos_mis = [brid for brid in brpos if brid != b38]
            brpos_hit = [brid for brid in brpos if brid == b38]
            ompos_mis = [omid for omid in ompos if omid != b38]
            ompos_hit = [omid for omid in ompos if omid == b38]
            if brpos_mis:
                if brpos_hit:
                    report.write('br map table has coords for %s that do not match with dbsnp %s (%s), but also coords that do (%s). No action coded yet\n' % (mid,b38,','.join(brpos_mis),','.join(brpos_hit)))
                else:
                    report.write('br map table coords for %s (%s) do not match dbsnp coord (%s). Correcting the coordinates in map table\n' % (mid,','.join(brpos_mis),b38))
                    [self.mtab_change_pos(mid,oldpos=mis,newpos=b38) for mis in brpos_mis]
            if ompos_mis:
                print('unchecked: omics b38 coord(s) for id %s (%s) do not match with the dbsnp coord (%s). Flagging entry in db\n' % (mid,','.join(ompos_mis),b38))
                [self.pos_flag(mid,chrpos,fl=-5) for chrpos in ompos_mis]
            if not ompos_hit:
                print('unchecked: dbsnp coord for id %s (%s) is not in omics db. Adding it.\n' % (mid,b38))
                self.addpos(mid,chrpos=b38,ds='dbsnp',build='38')
        report.close()
