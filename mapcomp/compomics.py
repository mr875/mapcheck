class CompOmics:

    def step(self,omics="omics",start=1,finish=None): # s&f: 1/0,5 -> 6,10 etc ...
        self.connectomics(omics)
        counter = 0
        try:
            if not finish:
                finish = self.qf.row_count
            for line in self.qf.readls():
                counter += 1
                if counter < start:
                    continue
                if counter == (finish + 1):
                    break
                print(line)
                self.tree(*self.getcols(line))
        finally:
            self.finish()

    def uid_proc(self,uid):
        print("'%s' found" % (uid))

    def tree(self,uid,suid,chrm,pos):
        q = 'SELECT EXISTS (SELECT * FROM consensus WHERE id = %s)'
        vals = (uid,)
        self.omcurs.execute(q,vals) 
        if self.omcurs.fetchone()[0]:
            self.uid_proc(uid)
            return
        if suid:
            vals = (suid,)
            self.omcurs.execute(q,vals) 
            if self.omcurs.fetchone()[0]:
                self.uid_proc(suid)
                return

