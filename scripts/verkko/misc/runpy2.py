# Richard Darst, November 2011
# Based on something from at least before December 2010

import runpy

def main(argv=None):
    """Simulate 'python -m <modname>'.

    This is like the ``runpy`` standard library module, but this
    allows chained -m invocations.

    argv: The argv of the program, default sys.argv

        argv

    Some notes:
    'python -m somemodule' =>
        argv[0]=the module calling this
    'python -m somemodule' arg =>
        argv[0]=module calling this
        argv[1]=arg
    """
    if argv is None:
        import sys
        argv = sys.argv

    print argv
    del argv[0] # The script calling this.

    if len(argv) > 0 and argv[0] == '-m':
        modname = argv[1]
        del argv[0:1]
        runpy.run_module(modname, run_name='__main__', alter_sys=True)
    elif len(argv) > 0:
        runpy.run_path(argv[0], run_name='__main__')

    else:
        from code import interact
        interact(local={'__name__':'__main__'})

if __name__ == "__main__":
    # Silly function.  This ignores each "-m runpy2"
    main()
