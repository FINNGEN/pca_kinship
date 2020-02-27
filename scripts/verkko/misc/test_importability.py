import importlib
import os
import re
import sys
import unittest

from nose.tools import *

"""
Tests the importability of all .py files in verkko.

It makes a list of every .py file in verkko, excluding files or
directories in the ``exclude_files`` variable.  It will then try to
import every one of them, with the DISPLAY environment variable
unset.
"""

# Files to exclude.  Prefix in "./".  Filenames are removed, directory
# names cause all included files to be excluded.
exclude_files = ['./docs', './misc/interactnow.py']

verkko_dir = os.path.relpath(os.path.dirname(os.path.dirname(__file__)))


def check_importability(fname):
    """This is the actual test function"""

    # Set sys.path
    old_syspath = sys.path
    sys.path.insert(0, '..')
    # remove DISPLAY from environment
    old_display = os.environ.get('DISPLAY', None)
    if 'DISPLAY' in os.environ:
        del os.environ['DISPLAY']

    # Get the true module name and do the actual test
    mname = fname[2:][:-3].replace('/', '.')
    #print ' ', mname
    importlib.import_module('verkko.'+mname)

    # Revert sys.path and DISPLAY to previous values
    sys.path = old_syspath
    if old_display is not None:
        os.environ['DISPLAY'] = old_display


def test_files():
    """Test generator, for all filenames.

    Calls check_importability on every filename."""
    for fname in all_python_files:
        yield check_importability, fname


# The function and os.path.walk call below create the list of all
# filenames.
all_python_files = [ ]
def callback(arg, dirname, fnames):
    for fname in fnames[:]:  # [:] to make a copy so we can modify in place
        fpath = os.path.join(dirname, fname)
        # Check exclusions
        if fpath in exclude_files:
            # removing from fnames means the directory won't be walked.
            fnames.remove(fname)
            continue
        # exclude editor backup files (emacs)
        if '#' in fname:
            continue
        # If it's a py file, add it to all_python_files.
        if fname.endswith('.py'):
            new_name = re.sub(r'[^a-zA-Z1-9_]+', r'_',
                              os.path.join(dirname, fname)).strip('_')

            all_python_files.append(fpath)
        #elif os.path.isdir(fname):
        #    if not os.path.exists(os.path.join(dirname,fname,'__init__.py')):
        #        fnames.remove(fname)

os.path.walk(verkko_dir, callback, None)
