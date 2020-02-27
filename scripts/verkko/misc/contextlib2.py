from __future__ import absolute_import

import contextlib
import os
import shutil
import tempfile


@contextlib.contextmanager
def chdir_context(dirname):
    """Context manager for chdir.

    Chdir to ``dirname``.  At exit, chdir to the former working
    directory.


    Arguments:
        dirname: str, directory to change to.


    Example::

        with chdir_context(self.dir):
            some_code
    """
    olddir = os.getcwd()
    os.chdir(dirname)
    yield
    os.chdir(olddir)

@contextlib.contextmanager
def tmpdir_context(chdir=False, delete=True, suffix='', prefix='tmp',dir=None):
    """Context manager for temporary directories.

    The context argument is the tmpdir name, absolute path.

    Arguments:

        chdir: bool, default False
            if true, chdir to the tmpdir when in the context.

        delete: bool, default True
            if true, delete the tmpdir on clean-up.

        suffix, prefix, dir: passed to ``tempfile.mkdtemp``.

    Usage::

        with tmpdir_context(dir='/tmp/') as dirname:
            f = open(os.path.join(dirname, 'tmp.txt'), 'w')
    """
    olddir = os.getcwd()
    tmpdir = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
    if chdir:
        os.chdir(tmpdir)
    yield tmpdir
    if chdir:
        os.chdir(olddir)
    if delete:
        shutil.rmtree(tmpdir)
