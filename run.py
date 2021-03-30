#from mapcomp.maptable import MapTable
from mapcomp.ce import CE
from mapcomp.ce import omni
from mapcomp.ce import ukb

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
#    msex = CE('msexome_map')
#    msex.add_file()
#    msex.step('chip_comp',0)
#    om = omni('omniexpress_map').table_dump()
#    om = omni('omniexpress_map')
#    om.add_file()
#    om.step('chip_comp',0)
#    om21 = omni('omniexpress_v2_1_map').table_dump()
#    om21 = omni('omniexpress_v2_1_map')
#    om21.add_file()
#    om21.step('chip_comp',0)
#    ukb21 = ukb('ukbbaffy_v2_1_map').table_dump()
#    ukb21 = ukb('ukbbaffy_v2_1_map')
#    ukb21.add_file()
#    ukb21.step('chip_comp',363835)
# should be possible now:
    ukb21 = ukb('ukbbaffy_v2_1_map','ukbb2_1_map_0121.txt')
    ukb21.add_file()
    ukb21.step('cc2',start=0,finish=None)

if __name__ == "__main__":
    main()

