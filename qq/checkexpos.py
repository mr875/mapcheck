# added qq/checkexpos.py to mapcheck code in case I want to set up some routines that look for correct positions for rsids in build 37. could think about inserting buuld 38 positions where they don't exist. Could also think about doing it for ALL rsids in db by splitting the into pieces and doing progressively. This is to deal with $proc/mr875/tasks/map/type_two inconsistencies (mostly build 37)
# to use this script source mapcheck/pth
import sys
from utils.connect import DBConnect

def main(argv):
    print(sys.path)

if __name__ == "__main__":
    main(sys.argv[1:])
