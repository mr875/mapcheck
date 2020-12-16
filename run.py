#from mapcomp.maptable import MapTable
from mapcomp.ce import CE

def main():

#    ce = CE('coreexome_map').table_dump()
    ce = CE('coreexome_map')
    ce.add_file() # if file is present already (to avoid running table_dump() every time during development)
    ce.step('cc4',0,10)

if __name__ == "__main__":
    main()

