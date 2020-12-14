import time

class CompOmics:

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
        now = int(time.time() - seconds)
        print("got to end, %s seconds" % (now))

    def uid_proc(self,uid):
        pass
        #print("'%s' found" % (uid))

    def tree(self,uid,suid,chrm,pos,alts):
        q = 'SELECT EXISTS (SELECT * FROM consensus WHERE id = %s)'
        vals = (uid,)
        self.omcurs.execute(q,vals) 
        if self.omcurs.fetchone()[0]:
            self.uid_proc(uid)
            return
        if alts:
            for altrs in alts:
                vals = (altrs,)
                self.omcurs.execute(q,vals) 
                if self.omcurs.fetchone()[0]:
                    self.uid_proc(altrs)
                    return
        if suid:
            vals = (suid,)
            self.omcurs.execute(q,vals) 
            if self.omcurs.fetchone()[0]:
                self.uid_proc(suid)
                return
        q = 'SELECT id FROM alt_ids WHERE alt_id = %s'
        vals = (uid,)
        self.omcurs.execute(q,vals) 
        result = [s[0] for s in self.omcurs.fetchall()]
        if len(result) > 0:
            if len(result) > 1:
                print("more than 1 possible id for ",uid,suid,chrm,pos,' : ',', '.join(result))
            self.uid_proc(result[0])
            return
        print("not found in consensus or alt_ids: ",uid,suid,chrm,pos)

