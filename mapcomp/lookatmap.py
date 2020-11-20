from connect import DBConnect

class LookMap:

    tabcols = {'coreexome_map':['snp','dbsnpid','chr']}

    def __init__(self,tabname):
        self.brconn = DBConnect("br")
        self.brcurs = self.brconn.getCursor()
        self.omconn = DBConnect("cc4")
        self.omcurs = self.omconn.getCursor()
        self.tabname = tabname

    def getrows(self,limit=None):
        q = 'SELECT '
        for col in self.tabcols[self.tabname]:
            q = q + col + ', '
        q = q[:-2] + ' FROM ' + self.tabname
        if limit:
            q = q + ' LIMIT ' + str(limit)
        print(q)

    def test_conn(self):
        qbr = 'SELECT * FROM ' + self.tabname + ' LIMIT 5'
        qom = 'SELECT * FROM flank LIMIT 3'
        try:
            self.brcurs.execute(qbr)
            self.omcurs.execute(qom)
            resultbr = self.brcurs.fetchall()
            resultom = self.omcurs.fetchall()
        except:
            self.finish()
            raise
        #print(resultbr)
        #print(resultom)

    def finish(self):
        self.brconn.close()
        self.omconn.close()
        
