from connect import DBConnect

class LookMap:

    def __init__(self):
        self.brconn = DBConnect("br")
        self.brcurs = self.brconn.getCursor()
        self.omconn = DBConnect("cc4")
        self.omcurs = self.omconn.getCursor()

    def test_conn(self):
        qbr = 'SELECT * FROM coreexome_map LIMIT 5'
        qom = 'SELECT * FROM flank LIMIT 3'
        try:
            self.brcurs.execute(qbr)
            self.omcurs.execute(qom)
            resultbr = self.brcurs.fetchall()
            resultom = self.omcurs.fetchall()
        except:
            self.finish()
            raise
        print(resultbr)
        print(resultom)

    def finish(self):
        self.brconn.close()
        self.omconn.close()
        
