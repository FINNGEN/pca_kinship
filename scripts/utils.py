import time,sys,os,mmap,gzip,subprocess,binascii,logging
import numpy as np
from tempfile import NamedTemporaryFile
from functools import partial
import multiprocessing,csv 
cpus = multiprocessing.cpu_count()

mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
mem_mib = mem_bytes/(1024.**2) 
proc_mem = mem_mib / (cpus +1)




log_levels = {
    'critical': logging.CRITICAL,
    'error': logging.ERROR,
    'warn': logging.WARNING,
    'warning': logging.WARNING,
    'info': logging.INFO,
    'debug': logging.DEBUG
}


def return_open_func(f):
    '''
    Detects file extension and return proper open_func
    '''
   

    file_path,file_root,file_extension = get_path_info(f)

    if 'bgz' in file_extension:
        #print('gzip.open with rb mode')
        open_func = partial(gzip.open, mode = 'rb')
    
    elif 'gz' in file_extension:
        #print('gzip.open with rt mode')
        open_func = partial(gzip.open, mode = 'rt')

    else:
        #print('regular open')
        open_func = open      
    return open_func

def progressBar(value, endvalue, bar_length=20):
    '''
    Writes progress bar, given value (eg.current row) and endvalue(eg. total number of rows)
    '''

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()



def identify_separator(f):
    open_func = return_open_func(f)
    with open_func(f) as i:header = i.readline().strip()
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(header)
    return dialect.delimiter
    
def timing_function(some_function):

    """
    Outputs the time a function takes  to execute.
    """

    def wrapper(*args,**kwargs):
        t1 = time.time()
        some_function(*args)
        t2 = time.time()
        print("Time it took to run the function: " + str((t2 - t1)))

    return wrapper

def get_filepaths(directory):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.


    
def get_path_info(path):
    file_path = os.path.dirname(path)
    basename = os.path.basename(path)
    file_root, file_extension = os.path.splitext(basename)
    return file_path,file_root,file_extension
    

def file_exists(fname):
    '''
    Function to pass to type in argparse
    '''
    if os.path.isfile(fname):
        return str(fname)
    else:
        print(fname + ' does not exist')
        sys.exit(1)
 

def pretty_print(string,l = 30):
    l = l-int(len(string)/2)
    print('-'*l + '> ' + string + ' <' + '-'*l)
    

def mapcount(filename):

    if not os.path.isfile(filename):
        raise ValueError("File doesn't exist")
    
    try:
        return count_lines(filename)
    except:
        return 0
    
def count_lines(filename):
    '''
    Counts line in file
    '''
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    while readline():
        lines += 1
    return lines

class SliceMaker(object):
    '''
    allows to pass around slices
    '''
    def __getitem__(self, item):
        return item

def line_iterator(f,separator=None,count = False,columns = SliceMaker()[:],dtype = str,skiprows = 0):
    '''
    Function that iterates through a file and returns each line as a list with separator being used to split.
    N.B. it requires that all elements are the same type
    '''

    if not separator:
        separator = identify_separator(f)

    open_func = return_open_func(f)
    i = open_func(f)
    for x in range(skiprows):next(i)

    if count is False:
        for line in i:
            yield np.array(line.strip().split(separator),dtype = str)[columns].astype(dtype)
    else:
        row = 0
        for line in i:
            row +=1 
            yield row,np.array(line.strip().split(separator),dtype = str)[columns].astype(dtype)



def basic_iterator(f,separator =None,skiprows = 0,count = False,columns = 'all'):
    '''
    Function that iterates through a file and returns each line as a list with separator being used to split.
    '''

    if not separator:
        separator = identify_separator(f)
        
    open_func = return_open_func(f)
    i = open_func(f)
    for x in range(skiprows):next(i)

    if count is False:
        for line in i:
            line =line.strip().split(separator)
            line = return_columns(line,columns)
            yield line
    else:
        row = 0
        for line in i:
            line =line.strip().split(separator)
            line = return_columns(line,columns)
            row += 1   
            yield row,line
def return_columns(l,columns):
    '''
    Returns all columns, or rather the elements, provided the columns
    '''
    if columns == 'all':
        return l
    elif type(columns) == int:
        return l[columns]
    elif type(columns) == list:
        return list(map(l.__getitem__,columns))

def tmp_bash(cmd,show=True,check = False,):

    if show:
        print(cmd)
    scriptFile = NamedTemporaryFile(delete=True)
    with open(scriptFile.name, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(cmd + "\n")

    os.chmod(scriptFile.name,0o777)
    scriptFile.file.close()

    if check:
        subprocess.check_call(scriptFile.name)
    else:
        subprocess.call(scriptFile.name,stderr = subprocess.DEVNULL)

def natural_sort(l):
    import re
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


def pad(s):
    '''
    Prepend/append an empty space to the input.
    '''
    return ' ' + str(s) + ' '

def return_header(f):
    open_func = return_open_func(f)
    with open_func(f) as i:header = i.readline().strip()
    delimiter = identify_separator(f)
    header = header.split(delimiter)
    return header
        
        
def valid_string(s):
    if s:
        return s
    else:
        print('Invalid String')
        sys.exit(1)
        

def merge_files(o_file,file_list):
    with open(o_file,'wt') as o:
        for f in file_list:
            with open(f,'rt') as i:
                for line in i:
                    o.write(line)

def make_sure_path_exists(paths):
    if type(paths) is not list: paths = [paths]
    import errno
    for path in paths:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise



def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'

def write_fam_samplelist(o_file,sample_list):
    with open(o_file,'wt') as o :
            for sample in sample_list:
                o.write(sample + '\t' + sample + '\n')
    
    

class Logger(object):
    "Lumberjack class - duplicates sys.stdout to a log file and it's okay"
    #source: https://stackoverflow.com/q/616645
    def __init__(self, filename="Red.Wood", mode="a"):
        self.stdout = sys.stdout
        self.file = open(filename, mode)
        sys.stdout = self

    def __del__(self):
        self.close()

    def __enter__(self):
        pass

    def __exit__(self, *args):
        self.close()

    def write(self, message):
        self.stdout.write(message)
        self.file.write(message)

    def flush(self):
        self.stdout.flush()
        self.file.flush()
        os.fsync(self.file.fileno())

    def close(self):
        if self.stdout != None:
            sys.stdout = self.stdout
            self.stdout = None

        if self.file != None:
            self.file.close()
            self.file = None



def plot_stacked_bar(data, series_labels, category_labels=None, 
                     show_values=False, value_format="{}", y_label=None, 
                     colors=None, grid=True, reverse=False):
    """Plots a stacked bar chart with the data and labels provided.

    Keyword arguments:
    data            -- 2-dimensional numpy array or nested list
                       containing data for each series in rows
    series_labels   -- list of series labels (these appear in
                       the legend)
    category_labels -- list of category labels (these appear
                       on the x-axis)
    show_values     -- If True then numeric value labels will 
                       be shown on each bar
    value_format    -- Format string for numeric value labels
                       (default is "{}")
    y_label         -- Label for y-axis (str)
    colors          -- List of color labels
    grid            -- If True display grid
    reverse         -- If True reverse the order that the
                       series are displayed (left-to-right
                       or right-to-left)
    """

    from matplotlib import pyplot as plt

    ny = len(data[0])
    ind = list(range(ny))

    axes = []
    cum_size = np.zeros(ny)

    data = np.array(data)

    if reverse:
        data = np.flip(data, axis=1)
        category_labels = reversed(category_labels)

    for i, row_data in enumerate(data):
        #axes.append(plt.bar(ind, row_data, bottom=cum_size, label=series_labels[i], color=colors[i]))
        axes.append(plt.bar(ind, row_data, bottom=cum_size, label=series_labels[i]))
        cum_size += row_data

    if category_labels:
        plt.xticks(ind, category_labels)

    if y_label:
        plt.ylabel(y_label)

    plt.legend()

    if grid:
        plt.grid()

    if show_values:
        for axis in axes:
            for bar in axis:
                w, h = bar.get_width(), bar.get_height()
                plt.text(bar.get_x() + w/2, bar.get_y() + h/2, 
                         value_format.format(h), ha="center", 
                         va="center")



  # sort inf type by kinship degree
    
    # labels = [k for k, v in sorted(degree_dict.items(), key=lambda item: item[1],reverse = True)]
    # print(labels)
    # for inftype in labels:
    #     #edges of specific type
    #     s = [(n,m) for (n,m) in G.edges if G[n][m]['weight'] == inftype]
    #     h = nx.Graph()
    #     h.add_edges_from(s)
    #     print(inftype,h.number_of_nodes(),h.number_of_edges())
    #     # degree distribution
    #     count = bins.bin_count([d for n,d in h.degree()])
    #     plot_data.append(count)
    
    # plot_data = np.array(plot_data)
    # # update legend labels
    # print(plot_data)
    # print(cat_labels)
    #plot_stacked_bar(plot_data,labels,category_labels= cat_labels)

