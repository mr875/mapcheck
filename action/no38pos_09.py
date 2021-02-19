import timeit
import re

# py fileproc.py out_files_corexome/out_sh/out_no38pos_rs.txt coreexome_map

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
        print('Done. %s seconds' % (int(timeit.default_timer() - self.grandstart)))
        report.close()
