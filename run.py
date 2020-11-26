from mapcomp.lookatmap import LookMap
from mapcomp.ce import CE

class LMCE(LookMap,CE):
    pass

def main():
    ce = LMCE('coreexome_map')
    ce.test_conn()
    try:
        ce.loadcurs_br(limit=5)
        ce.run_through()
    finally:
        ce.finish()

if __name__ == "__main__":
    main()

