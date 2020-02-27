import os.path

def get_dir_to_cur_module(file_attr):
    """
    Give the __file__ attribute of your current module
    and returns the absolute path to the directory in 
    which it is contained.

    Parameters
    ---------- 
    file_attr : str
        the __file__ attribute of the module you're executing

    Returns
    -------
    dir : str
        the aboslute path to the directory the caller module
        is in.
    """

    return os.path.dirname(os.path.abspath(file_attr))
