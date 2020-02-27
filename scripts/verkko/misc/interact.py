# David Creswick and Richard Darst, 2005

"""Improved code.interact()

Can take globals argument.  Can automatically find locals and globals
based on the stack - by default gets them from one stack frame up.
Suggested uses:

In Python, interact now::
  import fitz.interact ; interact.interact()

In Python, interact now.  This only works once, since the module is
only imported once::
  import fitz.interactnow

In Python, go into interact mode at a traceback::
  import fitz.interacttb

Run from command line and drop into PDB upon each traceback:
  python -m fitz.pdbtb
... or enable it from python:
  import fitz.pdbtb ; pdbtb.enable()
"""

import code
import inspect
import os
import sys

# Set up readline.  Perhaps this will interact with other things badly
# (i.e., if readline shouldn't be here), but I'll deal with that when
# it comes.
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")

scrollback = 20
line_continuers = (",", "\\")
tab = "        "

def interact(banner="", local=None, stackLevel=1, frame=None,
             mode='interact', history=False):
    """Interact using the current global and local scope and history.

    arguments:
    banner -- print this text before interacting.
    local -- ignored
    history: bool or int
        if non-false, this many lines of history are found from the
        file being interacted with and inserted into the readline
        buffer


    """
    if len(banner) > 0:
        print banner
    # Get data from calling frame
    #calling_data = inspect.stack()[stackLevel]
    #filename = calling_data[1]
    #lineno   = calling_data[2]
    #frame = calling_data[0]
    #local = frame.f_locals
    #global_ = frame.f_globals

    if frame is None:
        frame = inspect.stack()[stackLevel][0]
    filename = frame.f_code.co_filename
    lineno = frame.f_lineno
    local = frame.f_locals
    global_ = frame.f_globals


    if history:
        new_history = [ ]
        if os.access(filename, os.R_OK):
            new_history.extend(get_history(filename, lineno))
        # add lines to history, removing proper amount of whitespace
        for line in new_history:
            readline.add_history(line)

    # add an exit() function under the name __exit.  Hopefully this
    # won't cause namespace collisions!
    local["__exit"] = sys.exit
    if mode == 'pdb':
        import pdb
        p = pdb.Pdb()
        def do_quit(self, arg):
            sys.exit()
        p.do_quit = type(p.do_quit)(do_quit, p, pdb.Pdb)
        p.do_q = p.do_quit
        p.reset()
        p.interaction(frame, None)
        #from fitz import interactnow
        #pdb.post_mortem(frame)
    else:
        interact2(local=local, global_=global_, banner="")

def pdb(*args, **kwargs):
    return interact(*args, mode='pdb',
                    stackLevel=kwargs.pop('stackLevel', 1)+1,
                    **kwargs)

def get_history(filename, lineno):
    # get lines
    f = file(filename)
    lines = f.readlines()
    f.close()

    #get amount of leading whitespace
    current_line = lines[lineno-1]
    whitespace = len(current_line) - len(current_line.lstrip())

    # Get just the lines we want
    if lineno-scrollback >= 0:
        lines = lines[lineno-scrollback:lineno-1]
    else:
        lines = lines[0:lineno-1]
    # Reconstruct list of lines, taking into account wrapping
    newlines = [ ]
    curline = ""
    continuation = False
    for line in lines:
        line = line.replace("\t", tab)
        line = line.rstrip()
        if line == "":
            continue
        if continuation:
            line = line.lstrip()

        # Determine if we have a complete command.  If it is, append
        # it to newlines, and set continuation = False.  If it isn't,e
        # append it to curline and set continuation = True
        if len(line) > 0 and line[-1] == "\\":
            curline = "".join((curline, line[:-1]+" "))
            continuation = True
        elif len(line) > 0 and line[-1] == ",":
            curline = "".join((curline, line+" "))
            continuation = True
        else:
            newlines.append("".join((curline, line)))
            curline = ""
            continuation = False
    # Strip leading whitespace
    for i, line in enumerate(newlines):
        if len(line[:whitespace].strip()) == 0:
            # if it starts with just whitespace
            newlines[i] = line.rstrip()[whitespace:]
        else:
            # we need less whitespace here
            newlines[i] = line.strip()
    return newlines

def null_function(a, b):
    pass


# Hack up an interactive interperter. I just do the minimal
# modifications to what is in code.py.  Note, that I do use some
# Python distribution code here.  It is GPL compatible.

# The usage of the variable names `local` and `locals`, etc, is
# unfortunante.  We see both of them.  I try to be consistent with
# what already exists, but `global` is a statement!

class InteractiveConsole(code.InteractiveConsole):
    def __init__(self, locals=None,
                 globals=None,
                 filename="<console>"):
        """Documentation comes from `code.InteractiveConsole`"""
        # We use the previous initialization method, but now we also
        # set self.globals to what we passed.
        code.InteractiveConsole.__init__(self, locals)
        self.globals = globals

    def runcode(self, code_):
        # This comes straight from code.py
        try:
            #exec code in self.locals
            exec code_ in self.globals, self.locals
        except SystemExit:
            raise
        except:
            self.showtraceback()
        else:
            if code.softspace(sys.stdout, 0):
                print

def interact2(banner=None, readfunc=None, local=None, global_=None):
    """Closely emulate `code.interact`

    This function is exactly like `code.interact` (see it's
    documentation), but has the additional argument:
    globals -- mapping to use for global namespace namespace.
    """
    console = InteractiveConsole(locals=local,
                                 globals=global_)
    if readfunc is not None:
        console.raw_input = readfunc
    else:  # this is redundant since we imported readline above.
        try:
            import readline
        except ImportError:
            pass
    console.interact(banner)



if __name__ == "__main__":
    a = 12
#    print vars()
    if True:
        pass
    null_function(1, \
                  2)
    null_function(1, 2 +
                  2)
    null_function(1,
                  2)
    interact()

    print vars()
