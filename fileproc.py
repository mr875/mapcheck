import sys
import configparser
from utils.connect import DBConnect
from action.action import ProcFile
# py fileproc.py out_files_corexome/out_sh/out_map_new_alt_rs.txt coreexome_map

def loadConf():
    config = configparser.ConfigParser()
    config.read('conf/config.ini')
    return {'br':config['db_map']['br'],'omics':config['db_map']['omics']}

def main(argv):
    dbmap = loadConf()
    omics = DBConnect(dbmap['omics'])
    br = DBConnect(dbmap['br'])
    omcurs = omics.getCursor()
    brcurs = br.getCursor()
    try:
        fname = argv[0]
        tabname = argv[1]
    except IndexError:
        sys.exit("provide file_name and map table as arguments.eg:\npy fileproc.py out_files_corexome/out_sh/out_map_new_alt_rs.txt coreexome_map")
    repmode = True
    if len(argv) == 3:
        if 'False' in argv[2]:
            repmode = False # add switch to command line to execute db commands
    try:
        pf = ProcFile(fname,tabname,omcurs,brcurs,repmode)
    finally:
        omics.commit()
        omics.close()
        br.close()

if __name__ == "__main__":
    main(sys.argv[1:])
