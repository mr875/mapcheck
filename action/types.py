import re

class Types:

    maketable = 'CREATE TABLE IF NOT EXISTS extra_map (newid char(20) NOT NULL,linkid char(38) NOT NULL,chr char(4) DEFAULT NULL,GRCh38_pos int unsigned DEFAULT NULL,flank_seq varchar(1040) DEFAULT NULL,datasource varchar(40) NOT NULL,ds_chrpos varchar(20) DEFAULT NULL,chosen TINYINT DEFAULT 0,PRIMARY KEY (newid,linkid,datasource))'

    def newaltrs(self):
        #rs id not found in db but db already has another rs id for that variant
        report = open('report_newaltr_' + self.ts + '.txt',"w")
        self.extra_map_f = open('extra_map_' + self.ts + '.sql',"w")
        self.extra_map_f.write(self.maketable+';\n')
        count = 0
        for line in self.inp.read():
            count+=1
            if count == 60:
                break
            #print(line)
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
                self.swapout_main(swin=newrs,swout=currs)
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
                self.swapout_main(swin=newrs,swout=currs)
                continue
            newrsb38 = self.getb38(new_current[0]) 
            currsb38 = self.getb38(new_current[1])
            if newrsb38 != currsb38:
                ch_count = 0
                checkom = self.checkom_pos(mid=currs,chrpos=currsb38)
                chosen = 0
                ds_chrpos = None
                if not checkom[0]:
                    report.write('map table %s linked via alt id to omics db %s but they do not have the same chrpos in dbsnp. omics chrpos of %s does not match dbsnp. entry will be flagged in positions table. map table rs added to extra_map\n' % (newrs,currs,currs))
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
                if ch_count == 0 or ch_count == 2:
                   print('newrsb38 != currsb38, %s x reporting. check %s/%s' % (ch_count,newrs,currs)) 
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
                print('omics position for %s is wrong against dbsnp. dbsnp = %s. omics = %s. Entry to be flagged if not already' % (currs,currsb38,dbbd))
                self.pos_flag(mid=currs,chrpos=dbbd,fl=-5)
            if currsb38 not in dbpos:
                print('dbsnp position %s for omics id %s not in omics positions table. To be added' % (currsb38,currs))
                self.addpos(mid=currs,chrpos=currsb38,ds='dbsnp',build='38')
            for mapbp in mapbadpos:
                print('map table position for %s is wrong against dbsnp. dbsnp = %s. map table = %s.' % (newrs,currsb38,mapbp))
                self.mtab_change_pos(anid=newrs,oldpos=mapbp,newpos=currsb38)
            #if currsb38 not in dbpos: add it .... also newrs to be added as alt_id
            print("no action coded for map table %s and omics db %s" % (newrs,currs))
        report.close()
        self.extra_map_f.close()

