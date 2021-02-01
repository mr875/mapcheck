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
            print(mid,b38,merges)
        report.close()
