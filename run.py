from mapcomp.lookatmap import LookMap
from mapcomp.ce import CE

class LMCE(LookMap,CE):

    def run(self):
        self.connectomics('cc4')
        self.connectbr()
        self.test_conn()
        try:
            self.loadcurs_br(limit=5)
            self.run_through()
        finally:
            self.finish()

def main():

    ce = LMCE('coreexome_map').run()

if __name__ == "__main__":
    main()

