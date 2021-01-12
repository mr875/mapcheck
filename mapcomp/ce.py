from .maptable import MapTable
from .compomics import CompOmics

class CE(MapTable,CompOmics):
    
    tabcols = {'coreexome_map':['snp','dbsnpid','chr'],'humanexome_map':['snp','dbsnpid','chr'], 'infiniumimmunoarray_map':['snp','dbsnpid','chr'], 'msexome_map':['snp','dbsnpid','chr']}
    tabfilters = {'coreexome_map':['chr <> "chr"'], 'humanexome_map':['chr <> "chr"'],'msexome_map':['chr <> ""']}
    relvds = {'coreexome_map': ['170','171','232'],'humanexome_map':['233'],'infiniumimmunoarray_map':['70','234'],'msexome_map':['72']}

    def getcols(self,arr):
        uid,alts = self.rs_mult(arr[1]) #could be comma separated multiple, uid is expected in 2nd col (tabcols for CE type)
        suid = arr[0]
        if uid == '0' or uid == '':
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

class omni(MapTable,CompOmics):

    tabcols = {'omniexpress_map':['snp','dbsnpid','chr'],'omniexpress_v2_1_map':['snp','dbsnpid','chr']}
    tabfilters = {'omniexpress_map':['chr <> ""']}
    relvds = {'omniexpress_map':['73','235'],'omniexpress_v2_1_map':['73','235']}

    def getcols(self,arr):
        uid = arr[1]
        suid = arr[0]
        if "VG" in uid:
            uid = suid
            suid = None
        chrmpos = arr[2].split(':')
        if len(chrmpos) != 2:
            mssg = 'unexpected value in chr_pos field: ' + ', '.join(arr)
            raise RuntimeError(mssg)
        chrm = chrmpos[0]
        pos = chrmpos[1]
        return uid,suid,chrm,pos,[] #alts is an array

class ukb(CE):

    tabcols = {'ukbbaffy_v2_1_map':['chipid','dbsnpid','chr']}
    tabfilters = {'ukbbaffy_v2_1_map':['chr REGEXP "^[YX]|[0-9]+:"']}    
    relvds = {'ukbbaffy_v2_1_map':[231]}
