import re
# py fileproc.py out_files_corexome/out_sh/out_rsomics_rs.txt coreexome_map 

class OmRs_13:

    def omrs(self,brk=0,start=0):
        report = open('report_omrs_' + self.ts + '_' + self.tabname + '.txt',"w")
        count = 0
        for line in self.inp.read():
            count+=1
            if brk and count == brk:
                break
            if start and count < start:
                continue
            self.percent_comp(current=count,perc=10,total=brk,start=start)
            linesplit = re.split('\tlookup:',line.rstrip())
            mapid,omrs = re.split('\t',linesplit[0])
            mapid = re.split('=',mapid)[1]
            omrs = re.split('=',omrs)[1]
            stillmain = self.stillmain(omrs)
            if not stillmain:
                print('omics rs %s is not in consensus table any more, would have added it to map table id %s.\n' % (omrs,mapid))
                continue
            mgchk = self.mergecheck(linesplit[1])
            merges = mgchk['merges']
            becomes = self.mergedinto(evenarr=merges,potbef=omrs)
            mrgid = None
            if becomes: #omrs is merged to newer one
                mrgid = becomes[0]
                print('omics rs id %s is merged to %s and linked to map id %s\n' % (omrs,mrgid,mapid))
                mrgid_knwn = self.checkomid(cid=mrgid) # is the mrgid known? returns [conyesno,altyesno]
                if not mrgid_knwn[0] and not mrgid_knwn[1]:
                    print('merged id %s (from %s) is unknown to omics so it should be swapped in and added to map table id %s instead of %s\n' % (mrgid,omrs,mapid,omrs))
                else:
                    print('new merged id %s (from %s) is already known as consensus id or alternative id or both (%s). No action\n' % (mrgid,omrs,'-'.join([str(i) for i in mrgid_knwn])))
                    continue
            b38 = self.getb38(linesplit[1])
            if not b38:
                if not mrgid:
                    print('can not get dbsnp b38 position for omics rs id %s, map id %s. Skipping\n' % (omrs,mapid))
                else:
                    print('can not get dbsnp b38 position for new merge id %s (from %s). Adding as alternate id instead to omics, not adding to map table id %s' % (mrgid,omrs,mapid))
                    self.add_alt(alt=mrgid,main=omrs,ds='dbsnp')
                continue
            contig = False
            if '_' in b38: # can't edit chr:pos fields
                print('dbsnp coord is a contig reference for omics rs %s (or merge %s). Can not update position fields. map id %s\n' % (omrs,mrgid,mapid))
                contig = True
            getbrpos = self.checkbr_pos(mapid,b38)
            getompos = self.checkom_pos(omrs,chrpos=b38,build='38')
            if not getbrpos[0] and not contig:
                oldbrpos = [op for op in getbrpos[1] if op != '0:0' and op != 'pending' and op != 'flank_error']
                if len(oldbrpos) < 1:
                    report.write('map id %s, omics rs %s, possible merge %s. map table position not for changing/updating: %s\n' % (mapid,omrs,mrgid,','.join(getbrpos[1])))
                for op in oldbrpos:
                    print('map id %s, omics rs %s, possible merge %s. map position do not match dbsnp position (%s vs %s). Action: correct position in map table\n' % (mapid,omrs,mrgid,op,b38))
                    self.mtab_change_pos(anid=mapid,oldpos=op,newpos=b38)
            if not getompos[0] and not contig:
                print('omics rs id %s, possible merge %s, map table id %s, dbsnp position is unknown to omics, so adding it and flagging not matches (unless 0:0).\n' % (omrs,mrgid,mapid))                
                badom = [bp for bp in getompos[1] if bp != '0:0']
                self.addpos(omrs,chrpos=b38,ds='dbsnp',build='38')
                for bp in badom:
                    self.pos_flag(omrs,chrpos=bp,fl=-5)
            idin = omrs
            if mrgid:
                idin = mrgid
                print('omics rs %s due to be added to map table id %s but according to dbsnp it is merged to %s. So merged %s will be added to map table id %s and also swapped into omics db for %s\n' % (omrs,mapid,mrgid,mrgid,mapid,omrs))
                self.swapout_main(swin=mrgid,swout=omrs,ds='dbsnp')               
            else:
                print('omics rs %s to be added to map table id %s\n' % (omrs,mapid))
            self.mtab_change_id(xsting=mapid,chngto=idin) # NEED NEW METHOD TO PUT IN VS SWITCH
        report.close()
