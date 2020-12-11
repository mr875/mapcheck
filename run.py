#from mapcomp.maptable import MapTable
from mapcomp.ce import CE

def main():

    #ce = LMCE('coreexome_map').table_dump() 
    ce = CE('coreexome_map')
    ce.add_file() # if file is present already (to avoid running table_dump() every time during development)
    ce.step('cc4',0,5)

if __name__ == "__main__":
    main()

