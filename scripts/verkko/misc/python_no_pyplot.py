"""Wrapper to disable pylab/pyplot.

The indiscriminate use of pylab/pyplot is a major problem for
non-interactive use.  Commonly, this would be solved by running
``import matplotlib; matplotlib.use('Agg')``, but if you call this
multiple times (for example, in multiple test scripts), then it emits
a warning.

This module calls ``use('Agg')`` and then permanently deactivates the
``matplotlib.use`` function.

It is intented to be used as a script only, NOT as an importable
module.  Its only intended use so far is to

Instead of::

    python your_script.py

you can do::

    python -m verkko.misc.python_no_pyplot your_script.py

and it will work.  This uses the standard ``python -m`` mechanism.

"""

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    # Totally disable matplotlib.use
    matplotlib._use_orig = matplotlib.use
    matplotlib.use = lambda *args, **kwargs: None

    # Adjust sys.path and execute new script.
    import sys
    script = sys.argv[1]
    sys.argv = sys.argv[1:]

    # Compile code in namespace so that tracebacks are better.
    code = compile(open(script).read(), script, 'exec')
    exec(code)
