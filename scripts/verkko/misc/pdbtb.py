# Richard Darst, Nomember 2011

"""
"""
from __future__ import absolute_import


import inspect
import pdb
import sys

def excepthookpdb(t, value, tb):
    """Exception hook to invoke pdb on exceptions."""
    sys.__excepthook__(t, value, tb)
    pdb.post_mortem(tb)
def enable():
    """Enable exception hook."""
    sys.excepthook = excepthookpdb
def disable():
    """Disable exception hook."""
    sys.excepthook = sys.__excepthook__

def _run(func, *args, **kwargs):
    """Call .enable() and then run function."""
    enable()
    func(*args, **kwargs)

_run_history = [ ]
import readline
def _get_history(first=1):
    return [ readline.get_history_item(x)
             for x in range(first, readline.get_current_history_length()+1) ]
def _add_history(hist):
    [ readline.add_history(x) for x in hist ]
def _restore_history(hist):
    readline.clear_history()
    _add_history(hist)

def run(func, *args, **kwargs):
    """pdb hook: invokes pdb on exceptions in python.

    The function func is called, with arguments args and
    kwargs=kwargs.  If this func raises an exception, pdb is invoked
    on that frame.  Upon exit from pdb, return to python normally."""
    # save history
    old_hist = _get_history()
    old_hist_start = readline.get_current_history_length()+1

    try:
        return func(*args, **kwargs)
    except Exception as e:
        _add_history(_run_history)


        t, value, tb = sys.exc_info()
        sys.__excepthook__(t, value, tb)
        frame = sys.exc_info()[2]
        #tb = e.tb_frame
        pdb.post_mortem(tb)
        del frame   # the docs warn to avoid circular references.
        del t, value, tb

        _run_history[:] = _get_history(first=old_hist_start)
    readline.clear_history()
    _restore_history(old_hist)
    print old_hist
run_ipython = run

def now(frame=None, stackLevel=1):
    """Run pdb in the calling frame."""
    if frame is None:
        frame = inspect.stack()[stackLevel][0]
    p = pdb.Pdb()
    def do_quit(self, arg):
        # This raises SystemExit which escapes from pdb and returns to
        # normal execution.
        sys.exit()
    p.do_quit = type(p.do_quit)(do_quit, p, pdb.Pdb)
    p.do_q = p.do_quit
    p.reset()
    p.interaction(frame, None)


if __name__ == "__main__":
    from . import runpy2
    enable()
    runpy2.main()
