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

    def __init__(self, argv):
        dbmap = loadConf()
        omics = DBConnect(dbmap['omics'])
        self.omcurs = omics.getCursor()
        br = DBConnect(dbmap['br'])
        self.brcurs = br.getCursor()
        queryout = QueryFile('flankerr_out.txt',self.qry,db=dbmap['omics'])
        self.action_file(queryout)
        omics.close()
        br.close()

    def action_file(self,queryout):
        cnt = 0
        for mid in queryout.read():
            cnt += 1
            if cnt > 9:
                break

def main(argv):
    FlankErr(argv)

if __name__ == "__main__":
    main(sys.argv[1:])
