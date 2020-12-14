from .maptable import MapTable
from .compomics import CompOmics

class CE(MapTable,CompOmics):
    
    tabcols = {'coreexome_map':['snp','dbsnpid','chr']}
    tabfilters = {'coreexome_map':['chr <> "chr"']}
    relvds = {'coreexome_map': [170,171,232]}

    def getcols(self,arr):
        uid,alts = self.rs_mult(arr[1]) #could be comma separated multiple
        suid = arr[0]
        if uid == '0':
            uid = suid
            suid = None
        chrmpos = arr[2].split(':')
        if len(chrmpos) != 2:
            mssg = 'unexpected value in chr_pos field: ' + ', '.join(arr)
            raise RuntimeError(mssg)
        chrm = chrmpos[0]
        pos = chrmpos[1]
        return uid,suid,chrm,pos,alts #alts is an array

    def rs_mult(self,mval,dlim=','):
        arr = mval.split(dlim)
        if len(arr) < 2:
            return arr[0],[]
        nums = [int(rs.replace('rs','')) for rs in arr]
        smallest = min(nums) 
        ind = nums.index(smallest)
        chosen = arr.pop(ind)
        return chosen,arr

