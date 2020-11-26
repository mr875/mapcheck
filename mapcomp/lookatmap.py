import sys
from connect import DBConnect
#from .ce import CE

class LookMap:

    tabcols = {'coreexome_map':['snp','dbsnpid','chr']}

    def __init__(self,tabname):
        self.brconn = DBConnect("br")
        self.brcurs = self.brconn.getCursor()
        self.omconn = DBConnect("cc4")
        self.omcurs = self.omconn.getCursor()
        if tabname in self.related_names:
            self.tabname = tabname
        else:
            sys.exit("table name passed is '%s' but this name is not recognised for class '%s'" % (tabname,self.__class__.__name__))

    def loadcurs_br(self,limit=None):
        q = 'SELECT '
        for col in self.tabcols[self.tabname]:
            q = q + col + ', '
        q = q[:-2] + ' FROM ' + self.tabname
        if limit:
            q = q + ' LIMIT ' + str(limit)
        self.brcurs.execute(q)

    def run_through(self):
        for row in self.brcurs:
            print(row)

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
        
