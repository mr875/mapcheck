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
            new_current = re.split('\tlookup:',line.rstrip())
            mid,brtab_pos,brorig = re.split('\t',new_current[0])
            brorig = re.split('=',brorig)[1]
            if not self.stillmain:
                print('omics db id %s not present (in consensus table) any longer' % (mid))
                continue
            ompos = self.checkom_pos(mid=mid,build='38')[1] 
            if len(ompos) > 0:
                print('omics db id %s now does have b38 position. No action yet' % (mid))
                continue
            b38 = None
            becomes = []
            if len(new_current) > 1: # is dbsnp look up available 
                b38 = self.getb38(new_current[1])
                merges = self.mergecheck(new_current[1])['merges']
                if merges:
                    becomes = self.mergedinto(evenarr=merges,potbef=mid)
                print('dbsnp lookup available for %s. b38 = %s. merged into %s' % (mid,b38,','.join(becomes)))
            if not b38: # look up not available (and probably no rs id) OR b38 not available
                ompos37 = self.checkom_pos(mid=mid,build='37')[1] 
                print('omics id %s, map table coord %s, no dbsnp details available' % (mid,brtab_pos))
        print('Done. %s seconds' % (int(timeit.default_timer() - self.grandstart)))
        report.close()
