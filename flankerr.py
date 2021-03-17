import sys
import configparser
from utils.outfiles import OutFiles 
from utils.connect import DBConnect
from utils.queryfile import QueryFile, NormFile

brconn = None
brcurs = None
omconn = None
omcurs = None

def loadConf():                                                                                           
    config = configparser.ConfigParser()
    config.read('conf/config.ini')
    return {'br':config['db_map']['br'],'omics':config['db_map']['omics']}
            
def main(argv):
    dbmap = loadConf()
    omics = DBConnect(dbmap['omics'])
    qry = 'select distinct id from flank where id in (select id from flank group by id having sum(chosen) > 0) limit 10'
    queryout = QueryFile('flankerr_out.txt',qry,db=dbmap['omics'])

    omics.close()

if __name__ == "__main__":
        main(sys.argv[1:])
