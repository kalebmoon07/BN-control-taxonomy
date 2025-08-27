

def main():
    import bntaxonomy
    from bntaxonomy.iface import load_tools
    load_tools()
    # TODO
    return

if __name__ == "__main__":
    import sys
    import os.path
    libdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print(libdir)
    sys.path.insert(0, libdir)
    print(sys.path)
