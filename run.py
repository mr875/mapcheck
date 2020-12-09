from mapcomp.lookatmap import LookMap
from mapcomp.ce import CE

class LMCE(LookMap,CE):

    def run(self):
        try:
            self.connectomics('cc4')
            self.connectbr()
            self.test_conn()
            self.make_file(limit=5)
            self.run_through(limit=5)
        finally:
            self.finish()

def main():

    ce = LMCE('coreexome_map').run()

if __name__ == "__main__":
    main()

