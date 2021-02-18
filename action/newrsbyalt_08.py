import timeit
import re

# py fileproc.py out_files_corexome/out_sh/out_new_rs_byalt_rs.txt coreexome_map

class NewRsbyAlt_08:
# map table rsid not known to db, map table non-rsid in db alt_ids instead
    def newrsbyalt(self,brk=0,start=0):
        report = open('report_newrsbyalt_' + self.ts + '_' + self.tabname + '.txt',"w")
        count = 0
        for line in self.inp.read():
            count+=1
            if brk and count == brk:
                break
            if start and count < start:
                continue
            self.percent_comp(current=count,perc=10,total=brk,start=start)
            new_current = re.split('\tlookup= ',line.rstrip())
            mapid,linkid,midom = re.split('\t',new_current[0])
            mapid = re.split('= ',mapid)[1]
            linkid = re.split('= ',linkid)[1]
            midom = re.split('= ',midom)[1]
            mgchk = self.mergecheck(new_current[1])
            merges = mgchk['merges']
            print('map id = %s\tlink id = %s\tmain omics id = %s' % (mapid,linkid,midom))
            becomes = self.mergedinto(evenarr=merges,potbef=mapid)
            stillmain = self.stillmain(midom)
            mrgid = None
            if becomes:
                print('entered becomes. there\'s a merge to action')
                continue
                if len(becomes) > 1: # probably not necessary
                    print('map table id %s is merged into multiple (%s) according to dbsnp lookup. No action yet' % (mapid,'/'.join(becomes)))
                    continue
                mrgid = becomes[0]
                print('map table id %s merged to %s and linked to omics %s via %s. correcting map table (%s -> %s)\n' % (mapid,mrgid,midom,linkid,mapid,mrgid))
                self.mtab_change_id(xsting=mapid,chngto=mrgid)
                mrgid_knwn = self.checkomid(cid=mrgid) # is the new mrgid unknown as the mapid was?
                if not mrgid_knwn[0] and not mrgid_knwn[1]:
                    print('new id %s is also unknown to omics so adding old map id %s as alternative id to %s (if still present) and continuing\n' % (mrgid,mapid,midom))
                    if stillmain:
                        self.add_alt(alt=mapid,main=midom,ds=self.tabname)
                else:
                    print('new merged id %s is already known as consensus id or alternative id or both (%s). No action\n' % (mrgid,'-'.join([str(i) for i in mrgid_knwn])))
                    continue
            if not stillmain:
                print('main id %s (linked by %s) is not in consensus table anymore. Would have swapped it out with %s from map table (or merged %s)\n' % (midom,linkid,mapid,mrgid))
                continue
            b38 = self.getb38(new_current[1])
            if not b38:
                if not mrgid:
                    print('b38 position of map table id %s (omics %s, linked by %s) not available (may be withdrawn). Will add to omics as alternative id instead\n' % (mapid,midom,linkid))
                    self.add_alt(alt=mapid,main=midom,ds=self.tabname)
                else:
                    print('b38 position of new merge id %s (omics %s, linked by %s) not available (may be withdrawn). Will add to omics as alternative id instead\n' % (mrgid,midom,linkid))
                    self.add_alt(alt=mrgid,main=midom,ds='dbsnp')
                continue
            contig = False
            if '_' in b38:
                print('map table id %s (omics %s, linked by %s) dbsnp look up has a contig id instead of regular chromosome identifier: %s. Will swap %s (or merge %s) into omics for %s but will not correct/update positions\n' % (mapid,midom,linkid,b38,mapid,mrgid,linkid))
                contig = True

        print('Done. %s seconds' % (int(timeit.default_timer() - self.grandstart)))
