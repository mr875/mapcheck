import re
# py fileproc.py out_files_corexome/out_sh/out_rsomics_rs.txt coreexome_map 

class OmRs_13:

    def omrs(self,brk=0,start=0):
        report = open('report_omrs_' + self.ts + '_' + self.tabname + '.txt',"w")
        print("entered OmRs_13")
        count = 0
        for line in self.inp.read():
            count+=1
            if brk and count == brk:
                break
            if start and count < start:
                continue
            self.percent_comp(current=count,perc=10,total=brk,start=start)
            linesplit = re.split('\tlookup:',line.rstrip())
            mapid,omrs = re.split('\t',linesplit[0])
            mapid = re.split('=',mapid)[1]
            omrs = re.split('=',omrs)[1]
            mgchk = self.mergecheck(linesplit[1])
            merges = mgchk['merges']
            becomes = self.mergedinto(evenarr=merges,potbef=omrs)
            stillmain = self.stillmain(omrs)
            mrgid = None
            if becomes: #omrs is merged to newer one

        report.close()
