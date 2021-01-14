import sys
import configparser
from utils.connect import DBConnect
from action.action import ProcFile
# py fileproc.py out_files_corexome/out_sh/out_map_new_alt_rs.txt coreexome_map

def loadConf(self):
    config = configparser.ConfigParser()
    config.read('conf/config.ini')
    self.dbmap = {'br':config['db_map']['br'],'omics':config['db_map']['omics']}

def main(argv):
    try:
        fname = argv[0]
        tabname = argv[1]
    except IndexError:
        sys.exit("provide file_name and map table as arguments.eg:\npy fileproc.py out_files_corexome/out_sh/out_map_new_alt_rs.txt coreexome_map")
    pf = ProcFile(fname,tabname)

if __name__ == "__main__":
    main(sys.argv[1:])
