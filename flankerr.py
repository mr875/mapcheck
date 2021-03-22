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
    qvals = ()
    qfname = 'flankerr_out.txt'
    of = OutFiles('out_files')

    def __init__(self, argv):
        dbmap = loadConf()
        omics = DBConnect(dbmap['omics'])
        self.omcurs = omics.getCursor()
        br = DBConnect(dbmap['br'])
        self.brcurs = br.getCursor()
        queryout = QueryFile(self.qfname,self.qry,vals=self.qvals,db=dbmap['omics'])
        self.sqlout = self.of.new('out_' + self.__class__.__name__.lower() + '.sql')
        try:
            self.action_file(queryout)
        finally:
            self.of.closeall()
            omics.close()
            br.close()

    def get_dstomapt(self,mid):
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
        return mapts

    def action_file(self,queryout):
        cnt = 0
        for mid in queryout.read():
            cnt += 1
            mapts = self.get_dstomapt(mid)
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
            if not len(mapwhere) == len(mapval) == len(mapts_fltrd):
                raise Exception('lists mapwhere, mapval, mapts_fltrd should be the same size',mapwhere, mapval,mapts_fltrd)
            for ind,where in enumerate(mapwhere):
                self.printq(mapts_fltrd[ind],where,mapval[ind])
            #if cnt > 10:
             #   break
            
    def printq(self,maptable,where,vals):
        vals = ('flank_error',) + vals
        where = where.replace('%s','\'%s\'')
        full = 'UPDATE ' + maptable + ' SET chr = \'%s\' ' + where + ';\n'
        self.sqlout.write(full % vals)
        
class ChrPend(FlankErr):
    
    qry = 'SELECT id FROM positions where build = %s and chosen = %s and CONCAT(chr,\':\',pos) <> %s group by id,chosen having count(chosen) > 1'
    qvals = ('38',0,'0:0')
    qfname = 'chrpend_out.txt'

    def get_dstomapt(self,mid):
        allds = LinkID.get_positions_ds(mid,self.omcurs)
        nchose_ds = [var[0] for var in allds if var[1] == 0]
        mapts = []
        for ds in nchose_ds:
            if ds == '114':
                print('datasource DIL Taqman (114) for %s ignored' % (mid))
                continue
            if ds == 'dbsnp':
                logmess = 'mid %s has unresolved chr:pos despite a dbsnp entry (%s)'
                vals = (mid,','.join(nchose_ds))
                print(logmess % vals)
                self.sqlout.write('-- '+logmess % vals)
                self.sqlout.write('\n')
                continue
            mapts.extend(LinkID.anotds_to_mapt(ds))
        return mapts

    def printq(self,maptable,where,vals):
        vals = ('pending',) + vals
        where = where.replace('%s','\'%s\'')
        full = 'UPDATE ' + maptable + ' SET chr = \'%s\' ' + where + ';\n'
        self.sqlout.write(full % vals)

def main(argv):
    FlankErr(argv)
    #ChrPend(argv)


if __name__ == "__main__":
    main(sys.argv[1:])
