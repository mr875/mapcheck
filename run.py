#from mapcomp.maptable import MapTable
from mapcomp.ce import CE

def main():

#    ce = CE('coreexome_map').table_dump()
#    ce = CE('coreexome_map')
#    ce.add_file() # if file is present already (to avoid running table_dump() every time during development)
#    ce.step('chip_comp',0,20)
#    ce.step('chip_comp',0)
#    he = CE('humanexome_map').table_dump()
#    he = CE('humanexome_map')
#    he.add_file()
#    he.step('chip_comp',0,20)
#    he.step('chip_comp',0)
#    imuno = CE('infiniumimmunoarray_map').table_dump()
#    imuno = CE('infiniumimmunoarray_map')
#    imuno.add_file()
#    imuno.step('chip_comp',0)
#    msex = CE('msexome_map').table_dump()
    msex = CE('msexome_map')
    msex.add_file()
    msex.step('chip_comp',0)

if __name__ == "__main__":
    main()

