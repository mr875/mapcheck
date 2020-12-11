from .maptable import MapTable
from .compomics import CompOmics

class CE(MapTable,CompOmics):
    
    tabcols = {'coreexome_map':['snp','dbsnpid','chr']}
    relvds = {'coreexome_map': [170,171,232]}

    def getcols(self,arr):
        uid = arr[1]
        suid = arr[0]
        if uid == '0':
            uid = suid
            suid = None
        chrmpos = arr[2].split(':')
        if len(chrmpos) != 2:
            raise RuntimeError('unexpected value in chr_pos field %s' % (arr[2]))
        chrm = chrmpos[0]
        pos = chrmpos[1]
        return uid,suid,chrm,pos
