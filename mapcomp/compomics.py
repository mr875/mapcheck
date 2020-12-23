from utils.outfiles import OutFiles
from re import search
import time

class CompOmics:

    of = OutFiles('out_files')

    def output_setup(self):
        self.new_alt_rs = self.of.new(self.tabname + '_new_alt_rs.txt')
        self.new_rs_byalt = self.of.new(self.tabname + '_new_rs_byalt.txt')
        self.new_rs = self.of.new(self.tabname + '_new_rs.txt')
        self.rsomics = self.of.new(self.tabname + '_rs_from_omics.txt')
        self.OOpos = self.of.new(self.tabname + '_pos_for_oo.txt')
        self.new_mismatch = self.of.new(self.tabname + '_new_pos_mismatch.txt')
        self.no38pos = self.of.new(self.tabname + '_no38pos.txt')
        self.matchflag = self.of.new(self.tabname + '_matches_flagged.txt')
        self.matchflag_alt = self.of.new(self.tabname + '_matches_flagged_alt_avail.txt')
        self.lowcflank = self.of.new(self.tabname + '_lowconf_flank.txt')

    def step(self,omics="omics",start=1,finish=None): # s&f: 1/0,5 -> 6,10 etc ...
        self.output_setup()
        self.connectomics(omics)
        counter = 0
        seconds = int(time.time())
        try:
            if not finish:
                finish = self.qf.row_count
            todo = finish - start
            fiveper = int((todo/100)*5)
            fp_increment = fiveper
            fiveper = fiveper + start
            for line in self.qf.readls():
                counter += 1
                #print(counter,fiveper)
                if counter < start:
                    continue
                if counter == (finish + 1):
                    break
                if counter == fiveper:
                    done = int(((fiveper-start)/(todo))*100)
                    print('%s percent of %s (%s) lines parsed' % (done,todo,(fiveper-start)))
                    fiveper += fp_increment
                self.tree(*self.getcols(line))
        finally:
            self.finish()
            self.of.closeall()
        now = int(time.time() - seconds)
        print("got to end, %s seconds" % (now))

    def uid_fscan(self,uid,ori_uid):
        q = 'SELECT datasource,chosen,flank_seq FROM flank WHERE id = %s'
        vals = (uid,)
        self.omcurs.execute(q,vals)
        rows = self.omcurs.fetchall()
        dsrc = [row[0] for row in rows]
        chosen = [row[1] for row in rows]
        flank = [row[2] for row in rows]
        if 1 in chosen:
            for ind,ds in enumerate(dsrc):
                if chosen[ind] == 0 and ds in self.relvds[self.tabname]:
                    #print('low confidence flank detected\t%s (orig %s)\t%s' % (uid,ori_uid,flank[ind]))
                    self.lowcflank.write('low confidence flank detected\t%s (orig %s)\t%s\n' % (uid,ori_uid,flank[ind]))

    def uid_proc(self,uid,chrm,pos,ori_uid=None): # if ori_uid is used it means that uid is from omics db (via alt_ids table)
        if ori_uid:
            if self.isnewrs(uid,ori_uid):
                self.rsomics.write('mapfile id: %s\tomics db id: %s\n' % (ori_uid,uid))
        else:
            ori_uid = uid
        self.uid_fscan(uid,ori_uid)
        q = 'SELECT CONCAT(chr,":",pos),build,datasource,chosen FROM positions WHERE id = %s'
        vals = (uid,)
        self.omcurs.execute(q,vals)
        rows = self.omcurs.fetchall()
        chrpos = [row[0] for row in rows]
        build =  [row[1] for row in rows]
        ds =  [row[2] for row in rows]
        chosen =  [row[3] for row in rows]
        match = False # matched against normal db entry
        match_f = False # matched against flagged db entry
        mismatch = False # mismatch against normal db entry
        mismatch_f = False # mismatch against flagged db entry
        if '38' in build:
            for ind,cp in enumerate(chrpos):
                strrow = [str(c) for c in rows[ind]]
                if cp != (chrm + ':' + pos):
                    if build[ind] == '37':
                        continue
                    if chrm + ':' + pos == '0:0':
                        self.OOpos.write('positions available for %s 0:0\t%s\n' % (ori_uid,','.join(strrow)))
                        continue # SCENARIO 1
                    if chosen[ind] < 0:
                        mismatch_f = True
                        #print('mismatch against already flagged omics entry: %s (orig %s)\t%s:%s\t%s' % (uid,ori_uid,chrm,pos,','.join(strrow)))
                        continue # SCENARIO 2 (no report)
                    mismatch = True
                else: # match
                    if chosen[ind] < 0:
                        self.matchflag.write('%s (orig %s) matches flagged entry\t%s:%s\t%s\n' % (uid,ori_uid,chrm,pos,','.join(strrow)))
                        match_f = True # SCENARIO 3
                    else:
                        match = True # SCENARIO 4 (no report) 
            if mismatch:
                if not match_f: # full mismatch (against everything both flagged and not), 'a new mismatch' SCENARIO 5
                    self.new_mismatch.write('mismatch against all b38 db entries\t%s (orig %s)\t%s:%s\t%s\n' % (uid,ori_uid,chrm,pos,'/'.join(chrpos)))
                else: #SCENARIO 6
                    self.matchflag_alt.write('%s (orig %s) matched flagged, non flagged available\t%s:%s\t%s\n' % (uid,ori_uid,chrm,pos,','.join(chrpos )))
            if mismatch_f:
                if not match: #SCENARIO 7
                    print('no match, new to db\t%s (orig %s)\t%s:%s' % (uid,ori_uid,chrm,pos))
        else: #SCENARIO 8
            self.no38pos.write('no (b38) position available in omics:\t%s (orig %s)\t%s:%s\n' % (uid,ori_uid,chrm,pos))

    def fetchone(self,val,q):
        vals = (val,)
        self.omcurs.execute(q,vals) 
        return self.omcurs.fetchone()[0]

    def checkalt(self,val,q=None):
        if not q:
            q = 'SELECT id FROM alt_ids WHERE alt_id = %s'
        vals = (val,)
        self.omcurs.execute(q,vals) 
        result = [s[0] for s in self.omcurs.fetchall()]
        if len(result) > 0:
            result = list(dict.fromkeys(result))
            if len(result) > 1:
                print("more than 1 possible main id for ",val,' : ',', '.join(result))
            return result[0]
        return None

    def isnewrs(self,pid,sid):
        pid = search('^rs',pid)
        sid = search('^rs',sid)
        if pid and not sid: 
            return True
        return False

    def tree(self,uid,suid,chrm,pos,alts):
        q = 'SELECT EXISTS (SELECT * FROM consensus WHERE id = %s)'
        if self.fetchone(uid,q):
            self.uid_proc(uid,chrm,pos)
            return
        if alts:
            for altrs in alts:
                if self.fetchone(altrs,q):
                    self.uid_proc(altrs,chrm,pos)
                    return
        if suid:
            if self.fetchone(suid,q):
                if self.isnewrs(uid,suid): 
                    self.new_rs.write("not in consensus (uid): %s\tin consensus (suid): %s\n" % (uid,suid)) 
                self.uid_proc(suid,chrm,pos)
                return
        diffmain = self.checkalt(uid)
        if diffmain:
            self.uid_proc(diffmain,chrm,pos,uid)
            return
        if alts:
            for altrs in alts:
                diffmain = self.checkalt(altrs)
                if diffmain:
                    self.uid_proc(diffmain,chrm,pos,uid)
                    return
        diffmain = self.checkalt(suid)
        if diffmain:
            if self.isnewrs(uid,suid):
                if search('^rs',diffmain):
                    self.new_alt_rs.write('not in alt_ids (uid): %s\tin alt_ids (suid): %s\tdb main id: %s\n' % (uid,suid,diffmain))
                else:
                    self.new_rs_byalt.write('not in alt_ids (uid): %s\tin alt_ids (suid): %s\tdb main id: %s\n' % (uid,suid,diffmain))
            self.uid_proc(diffmain,chrm,pos,uid)
            return
        print("WARNING: not found in consensus or alt_ids: ",uid,suid,chrm,pos)

