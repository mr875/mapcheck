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
                self.tree(*self.getcols(line))
        finally:
            self.finish()

    def tree(self,uid,suid,chrm,pos):
        print(uid)
