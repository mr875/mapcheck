import sys
import configparser
from utils.outfiles import OutFiles 
from utils.connect import DBConnect
from utils.queryfile import QueryFile, NormFile
from action.linkid import LinkID

def loadConf():
    config = configparser.ConfigParser()
    config.read('conf/config.ini')
    return {'br':config['db_map']['br'],'omics':config['db_map']['omics']}
            
class FlankErr: 
    
    qry = 'select distinct id from flank where id in (select id from flank group by id having sum(chosen) > 0)'
    qfname = 'flankerr_out.txt'

    def __init__(self, argv):
        dbmap = loadConf()
        omics = DBConnect(dbmap['omics'])
        self.omcurs = omics.getCursor()
        br = DBConnect(dbmap['br'])
        self.brcurs = br.getCursor()
        queryout = QueryFile(self.qfname,self.qry,db=dbmap['omics'])
        try:
            self.action_file(queryout)
        finally:
            omics.close()
            br.close()

    def action_file(self,queryout):
        cnt = 0
        for mid in queryout.read():
            cnt += 1
            allds = LinkID.get_flank_ds(mid,self.omcurs)
            nchose_ds = [var[0] for var in allds if var[1] == 0]
            mapts = []
            for ds in nchose_ds:
                if ds == '114':
                    print('datasource DIL Taqman (114) for %s ignored' % (mid))
                    continue
                if ds == 'ext':
                    raise Exception('datasource %s for %s should not have a flank error\n' % (ds,mid))
                mapts.extend(LinkID.anotds_to_mapt(ds))
            #print(mapts)
            mapwhere = []
            mapval = []
            mapts_fltrd = []
            for mt in mapts:
                where,val,colwithid,res = LinkID.get_maptab_where(mid,mt,self.brcurs,self.omcurs)
                if res:
                    mapts_fltrd.append(mt)
                    mapwhere.append(where)
                    mapval.append(val)
            #print(mapts_fltrd,mapwhere,mapval)
            if not len(mapwhere) == len(mapval) == len(mapts_fltrd):
                raise Exception('lists mapwhere, mapval, mapts_fltrd should be the same size',mapwhere, mapval,mapts_fltrd)
            for ind,where in enumerate(mapwhere):
                full = 'UPDATE ' + mapts_fltrd[ind] + ' SET ... ' + where 
                print(full % mapval[ind])
            if cnt > 3:
                break

def main(argv):
    FlankErr(argv)

if __name__ == "__main__":
    main(sys.argv[1:])
