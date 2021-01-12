import sys
from action.action import ProcFile

def main(argv):
    try:
        fname = argv[0]
    except IndexError:
        sys.exit("provide file name as argument")
    pf = ProcFile(fname)

if __name__ == "__main__":
    main(sys.argv[1:])
