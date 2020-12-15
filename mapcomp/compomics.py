from utils.outfiles import OutFiles
from re import search
import time

class CompOmics:

    of = OutFiles('out_files')

    def step(self,omics="omics",start=1,finish=None): # s&f: 1/0,5 -> 6,10 etc ...
        self.connectomics(omics)
        counter = 0
        seconds = int(time.time())
        try:
            if not finish:
                finish = self.qf.row_count
            todo = finish - start
            fiveper = int((todo/100)*5)
            fp_increment = fiveper
            for line in self.qf.readls():
                counter += 1
                if counter < start:
                    continue
                if counter == (finish + 1):
                    break
                if counter == fiveper:
                    done = int((fiveper/todo)*100)
                    print('%s percent of %s (%s) lines parsed' % (done,todo,fiveper))
                    fiveper += fp_increment
                self.tree(*self.getcols(line))
        finally:
            self.finish()
            self.of.closeall()
        now = int(time.time() - seconds)
        print("got to end, %s seconds" % (now))

    def uid_proc(self,uid):
        pass
        #print("'%s' found" % (uid))

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
            self.uid_proc(uid)
            return
        if alts:
            for altrs in alts:
                if self.fetchone(altrs,q):
                    self.uid_proc(altrs)
                    return
        if suid:
            if self.fetchone(suid,q):
                if self.isnewrs(uid,suid): 
                    print("not in consensus (uid): %s\tin consensus (suid): %s" % (uid,suid))
                self.uid_proc(suid)
                return
        diffmain = self.checkalt(uid)
        if diffmain:
            self.uid_proc(diffmain)
            return
        if alts:
            for altrs in alts:
                diffmain = self.checkalt(altrs)
                if diffmain:
                    self.uid_proc(diffmain)
                    return
        diffmain = self.checkalt(suid)
        if diffmain:
            if self.isnewrs(uid,suid):
                print('not in alt_ids (uid): %s\tin alt_ids (suid): %s\tdb main id: %s' % (uid,suid,diffmain))
            self.uid_proc(diffmain)
            return
        print("WARNING: not found in consensus or alt_ids: ",uid,suid,chrm,pos)

