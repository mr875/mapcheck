import re

class Types:

    def percent_comp(self,current,perc,total=0):
        if not total:
            total = self.inp.row_count
        steps = int((perc/100.0)*total)
        #print('%s & %s' % (current,steps))
        try:
            if current % steps == 0:
                print('%s of %s (%s percent) done' % (current,total,int((current/total)*100)))
        except ZeroDivisionError:
            pass

    def newaltrs(self,brk=0):
        #rs id not found in db but db already has another rs id for that variant
        report = open('report_newaltr_' + self.ts + '.txt',"w")
        count = 0
        for line in self.inp.read():
            count+=1
            if brk and count == brk:
                break
            #print(line)
            self.percent_comp(current=count,perc=10,total=brk)
            new_current = re.split('Current',line)
            newrs = self.grabrsinp(new_current[0])
            currs = self.grabrsinp(new_current[1])
            newmchk = self.mergecheck(new_current[0])
            curmchk = self.mergecheck(new_current[1])
            if newmchk['withdrawn']:
                if not curmchk['withdrawn']:
                    report.write('map table %s is withdrawn from dbsnp, use omics db %s instead\n' % (newrs,currs))
                    self.add_alt(alt=newrs,main=currs,ds=self.tabname)
                    self.mtab_change_id(xsting=newrs,chngto=currs)
                else:
                    report.write('map table %s AND omics db %s are withdrawn from dbsnp\n' % (newrs,currs))
                continue   
            if curmchk['withdrawn']:
                report.write('omics db %s is withdrawn. map table %s to be used in omics db\n' % (currs,newrs))
                self.swapout_main(swin=newrs,swout=currs,ds=self.tabname)
                continue
            mergelist = newmchk['merges'] + curmchk['merges']
            newtocurr = self.frominto(evenarr=mergelist,bfor=newrs,aftr=currs)
            if newtocurr['correct']:
                report.write("map table %s merged into omics db %s, use %s\n" % (newrs,currs,currs))
                self.add_alt(alt=newrs,main=currs,ds=self.tabname)
                self.mtab_change_id(xsting=newrs,chngto=currs)
                continue
            currtonew = self.frominto(evenarr=mergelist,bfor=currs,aftr=newrs)
            if currtonew['correct']:
                report.write("omics db %s merged into map table %s, omics db to take %s\n" % (currs,newrs,newrs))
                self.swapout_main(swin=newrs,swout=currs,ds=self.tabname)
                continue
            newrsb38 = self.getb38(new_current[0]) 
            currsb38 = self.getb38(new_current[1])
            if not newrsb38 and not currsb38: #can't get dbsnp pos for both
                raise Exception('can not get dbsnp b38 coord for map table %s OR omics db %s. Uncoded scenario at the moment' % (newrs,currs))
            if not newrsb38 or not currsb38: #can't get dbsnp for one
                br_pos = self.checkbr_pos(newrs)[1]
                om_pos = self.checkom_pos(currs)[1]
                if not newrsb38: #have dbsnp pos for omics version
                    print(br_pos,om_pos)
                    bad_om_pos = [bp for bp in om_pos if bp != currsb38]
                    if len(bad_om_pos) == len(om_pos): # same as- if currsb38 not in om_pos:
                        print('unchecked innoncinnill: dbsnp position %s for omics id %s not in omics positions table. To be added\n' % (currsb38,currs))
                        self.addpos(mid=currs,chrpos=currsb38,ds='dbsnp',build='38')
                    for bp in bad_om_pos:
                        print('unchecked innoncinnfbib: omics position for %s is wrong against dbsnp. dbsnp = %s. omics = %s. Entry to be flagged if not already\n' % (currs,currsb38,bp))
                        self.pos_flag(mid=currs,chrpos=bp,fl=-5)
                    if currsb38 not in br_pos:
                        print('unchecked innoncificnib: position for map table id %s (%s) not found in dbsnp. Corresponding omics id %s, has its position available in dbsnp, as non matching coord %s. May not be the same variant. Adding to extra_map\n' % (newrs,','.join(br_pos),currs,currsb38))
                        self.extra_map(newid=newrs,linkid=currs,chrpos=None,datasource=self.tabname,chosen=-2,ds_chrpos=','.join(br_pos))
                    else:
                        report.write('position for map table id %s is not found in dbsnp. Corresponding omics id %s, has its position available in dbsnp which matches map table version. omics db id dbsnp pos = %s. map table id position = %s. Adding %s as alt_id\n' % (newrs,currs,currsb38,','.join(br_pos),newrs))
                        self.add_alt(alt=newrs,main=currs,ds=self.tabname)
                else: # not currsb38
                    pass
                continue
            if newrsb38 != currsb38:
                ch_count = 0
                checkom = self.checkom_pos(mid=currs,chrpos=currsb38)
                chosen = 0
                ds_chrpos = None
                if not checkom[0]:
                    report.write('map table %s linked via alt id to omics db %s but they do not have the same chrpos in dbsnp. omics chrpos of %s does not match dbsnp. entry will be flagged in positions table. dbsnp position %s will be added to positions table under %s. map table rs added to extra_map\n' % (newrs,currs,currs,currsb38,currs))
                    ch_count+=1
                    om_chrpos = list(dict.fromkeys(checkom[1]))
                    [self.pos_flag(currs,om_cp) for om_cp in om_chrpos]
                    self.addpos(mid=currs,chrpos=currsb38,ds='dbsnp')
                checkbr = self.checkbr_pos(mid=newrs,chrpos=newrsb38)
                if not checkbr[0]:
                    report.write('map table %s linked via alt id to omics db %s but they do not have the same chrpos in dbsnp. map table chrpos of %s does not match dbsnp. map table rs added to extra_map with error flag\n' % (newrs,currs,newrs))
                    ch_count+=1
                    chosen = -1
                    ds_chrpos = list(dict.fromkeys(checkbr[1]))
                    ds_chrpos=','.join(ds_chrpos)
                self.extra_map(newid=newrs,linkid=currs,chrpos=newrsb38,datasource=self.tabname,chosen=chosen,ds_chrpos=ds_chrpos)
                continue
            # if reached here then newrsb38 == currsb38 but they weren't found in either merge list so they are probably different at the variant type level
            dbpos = self.checkom_pos(mid=currs)[1]
            mapos = self.checkbr_pos(mid=newrs)[1]
            allpos = [currsb38] + dbpos + mapos
            allpos = list(dict.fromkeys(allpos))
            if len(allpos) == 1: # all chr:pos agree
                self.add_alt(alt=newrs,main=currs,ds=self.tabname)
                report.write('map table %s and omics %s have the same position on dbsnp and with each other. map table %s to be added as alt id to omics %s. They are likely to be a different variant-type\n' % (newrs,currs,newrs,currs))
                continue
            badpos = [ps for ps in allpos if ps != currsb38]
            dbadpos = [ps for ps in badpos if ps in dbpos]
            mapbadpos = [ps for ps in badpos if ps in mapos]
            for dbbd in dbadpos:
                report.write('omics position for %s is wrong against dbsnp. dbsnp = %s. omics = %s. Entry to be flagged if not already\n' % (currs,currsb38,dbbd))
                self.pos_flag(mid=currs,chrpos=dbbd,fl=-5)
            if currsb38 not in dbpos:
                report.write('dbsnp position %s for omics id %s not in omics positions table. To be added\n' % (currsb38,currs))
                self.addpos(mid=currs,chrpos=currsb38,ds='dbsnp',build='38')
            for mapbp in mapbadpos:
                report.write('map table position for %s is wrong against dbsnp. dbsnp = %s. map table = %s. map table coord to be changed/added\n' % (newrs,currsb38,mapbp))
                self.mtab_change_pos(anid=newrs,oldpos=mapbp,newpos=currsb38)
            self.add_alt(alt=newrs,main=currs,ds=self.tabname)
            report.write("omics id %s and map table id %s have the same coord in dbsnp. map table id %s to be added to omics id %s as alt_id\n" % (currs,newrs,newrs,currs))
        report.close()

