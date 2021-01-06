import sys
from utils.connect import DBConnect
from utils.queryfile import QueryFile, NormFile

class MapTable:

    brconn = None
    brcurs = None
    omconn = None
    omcurs = None
    qf = None #for QueryFile instance

    def __init__(self,tabname):
        if tabname in self.tabcols:
            self.tabname = tabname
            self.fname = 'out_' + self.tabname + '.txt'
        else:
            sys.exit("table name passed is '%s' but this name is not recognised for class '%s'" % (tabname,self.__class__.__name__))

    def table_dump(self,br="br",limit=None):
        self.make_file(br,limit)
        print("length of %s: %s" % (self.fname,self.qf.row_count))

    def connectbr(self,br="br"):
        if self.brconn:
            print("won't run %s.connectbr() because an existing connection already detected" % (self.__class__.__name__))
            return
        self.brconn = DBConnect(br)
        self.brcurs = self.brconn.getCursor()

    def connectomics(self,omics="omics"):
        if self.omconn:
            print("won't run %s.connectomics() because an existing connection already detected" % (self.__class__.__name__))
            return
        self.omconn = DBConnect(omics)
        self.omcurs = self.omconn.getCursor()

    def make_q(self,limit=None):
        q = 'SELECT DISTINCT '
        for col in self.tabcols[self.tabname]:
            q = q + col + ', '
        q = q[:-2] + ' FROM ' + self.tabname
        if self.tabname in self.tabfilters:
            flist = self.tabfilters[self.tabname] 
            fstr = ' WHERE ' + ' AND '.join(flist)
            q = q + fstr
        if limit:
            q = q + ' LIMIT ' + str(limit)
        return q

    def run_through(self,limit=None):
        self.brcurs.execute(self.make_q(limit))
        for row in self.brcurs:
            print(row)

    def make_file(self,dbname="br",limit=None):
        q = self.make_q(limit)
        self.qf = QueryFile(self.fname,q,(),dbname)

    def add_file(self):
        if self.qf:
            print("won't run %s.add_file() because a file is already in play" % (self.__class__.__name__))
            return
        self.qf = NormFile(self.fname)
        print("length of %s: %s" % (self.fname,self.qf.row_count))

    def test_conn(self):
        if not self.brcurs or not self.omcurs:
            print("activate db connections before running %s.test_conn()" % (self.__class__.__name__))
            return
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
        if self.brconn:
            self.brconn.close()
        if self.omconn:
            self.omconn.close()
        

