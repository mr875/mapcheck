import os
import datetime

class OutFiles:
    
    def __init__(self,directory=None):
        now = datetime.datetime.now().strftime("%b%d")
        self.path = os.getcwd()
        self.directory = directory
        self.fileobs = []
        if directory:
            self.path = self.path + '/' + directory + '_' + now
            if not os.path.exists(self.path):
                os.makedirs(self.path)

    def new(self,name):
        f = open(name,'w') # usage: f.write('text')
        self.fileobs.append(f)
        return f

    def closeall(self):
        for f in self.fileobs:
            f.close()

