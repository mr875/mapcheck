from mapcomp.lookatmap import LookMap

def main():
    gen = LookMap('coreexome_map')
    gen.test_conn()
    try:
        gen.getrows(5)
    finally:
        gen.finish()

if __name__ == "__main__":
    main()

