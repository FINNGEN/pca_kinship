from utils import basic_iterator,mapcount,np,subprocess,return_open_func,identify_separator,return_header,make_sure_path_exists,write_fam_samplelist
from collections import defaultdict as dd
from collections import Counter
import pandas as pd
import networkx as nx
import os.path,shlex,shutil

kinship_range = {0:[0.354,1],1: [0.177, 0.354], 2: [0.0884, 0.177], 3: [0.0442, 0.0884]}

def kinship(args):
    '''
    Builds a new plink file for king and performs king analysis
    '''
  
    # CALCULATE KINSHIP WITHIN DEGREE 3 if file is not passed
    if not args.kin:
        args.kin = args.kinPath + args.name +'.kin0'
        if not os.path.isfile(args.kin) or args.force:
            download_king()
            print('generating kinship relations file to ',args.kin)
            cmd = f'king -b {args.bed}  --kinship --degree 3 --cpus {args.cpus} --prefix {args.kin.split(".kin0")[0]}'
            subprocess.call(shlex.split(cmd))
        
    else:
        args.v_print(3,'kinship already calculated')

 
    # Return LIST OF RELATED INDIVIDUALS TO BE REMOVED AND LIST OF DUPLICATES
    args.related_couples = os.path.join(args.kinPath,args.name  +'_related_couples_' + str(args.degree) + '.txt')
    args.duplicates = os.path.join(args.kinPath,args.name  +'_duplicates.txt')
    if not os.path.isfile(args.related_couples) or not os.path.isfile(args.duplicates) or args.force:
        args.force = True
        print('Generating degree ' + str(args.degree) + ' related couples ... ')
        
        columns = [return_header(args.kin).index(elem) for elem in ['ID1','ID2','Kinship']]
        iterator = basic_iterator(args.kin,columns = columns,skiprows=1)
        
        duplicates = []
        deg_dict = dd(lambda : 3)
        with open(args.related_couples,'wt') as o:
            for sample1,sample2,kinship_value in iterator:
              
                kinship_deg = return_degree(float(kinship_value))
                if kinship_deg <= args.degree:
                    o.write(sample1 +'\t'+ sample2 + '\n')
                if kinship_deg ==0:
                    duplicates.append(sample1)
                    duplicates.append(sample2)
                for sample in sample1,sample2:
                    old_deg  = deg_dict[sample]
                    deg_dict[sample] = min(deg_dict[sample],kinship_deg)
            
                    
        duplicates = set(duplicates)
        write_fam_samplelist(args.duplicates,duplicates)
        print('done.')
        with open(os.path.join(args.kinPath,args.name  +'_kinship_count.txt'),'wt') as o :
            count = Counter(deg_dict.values())
            for i in range(0,4):
                o.write(str(i) + '\t' + str(count[i]) + '\n')
            
    else:
        pass
    
    print(f'degree {args.degree} related couples : {mapcount(args.related_couples)}')
    print(f"total duplicates : {mapcount(args.duplicates)}")

    # RETURN SET OF RELATED INDIVIDUALS WITH GREEDY AND NX ALGORITHMS
    args.related_individuals = args.kinPath + args.name + '_related_individuals_' + str(args.degree) + '.txt'
    if not  os.path.isfile(args.related_individuals) or args.force:
        args.force = True

        #import graph
        print(f'Generating degree {args.degree} related individuals ... ')
        print('Building nx graph..')
        g = nx.read_edgelist(args.related_couples)
        print(g.number_of_nodes())
        if args.test:
            samples = np.loadtxt(args.sample_fam ,usecols = [0],dtype = str)
            g = nx.Graph(g.subgraph(samples))

        print(g.number_of_nodes())
        print('Done.')
        
        # remove false finns
        print(f'removing non finns before starting to trim nodes {args.false_finns}')
        remove_nodelist = np.loadtxt(args.false_finns,dtype = str,usecols =[0])
        g.remove_nodes_from(remove_nodelist)
    
        # native nx algorithm vs greedy algorithm
        unrelated_nodes = []
        related_greedy = []
        for subgraph in connected_component_subgraphs(g):
            unrelated_nodes += nx.maximal_independent_set(subgraph)
            related_greedy += greedy_algorithm(subgraph)
            
        #sanity checks
        sanity_check(g,unrelated_nodes)
        sanity_check(g,g.nodes() - related_greedy)

        related_nodes = list(g.nodes() - set(unrelated_nodes))
        print(f'{len(related_nodes)} native related')
        print(f'{len(related_greedy)} greedy related')

        # return smaller set of related nodes from the two algorithms
        to_be_removed = min(related_greedy,related_nodes)             
        write_fam_samplelist(args.related_individuals,to_be_removed)
    else:
        args.v_print(1,'list of related samples already generated.')
        
    print(f'degree {args.degree} related individuals : {mapcount(args.related_individuals)}')
    

def sanity_check(graph,nodes):
    '''
    Given a list of nodes it makes sure that the algorithms are working properly.
    That is that the subgraph induced by the remainign nodes does not contain edges.
    '''

    assert graph.subgraph(nodes).number_of_edges() == 0

    
def greedy_algorithm(g):
    """
    Removes sequentially node with highest degree until there are no nodes left
    """
    
    degrees = dict(g.degree())
    removedNodes = []
    starting_edges = g.number_of_edges()

    edges = g.number_of_edges() 
    while edges >0:
        #find highest degree node
        maxNode = max(degrees, key=degrees.get)
        for neighbor in g[maxNode]:
            degrees[neighbor] -= 1
        removedNodes.append(maxNode)
        
        #delete node from degree dict and from network
        del degrees[maxNode]
        g.remove_node(maxNode)
        # update number of edges
        edges =g.number_of_edges()
        #progressBar(starting_edges - edges, starting_edges, bar_length=20)
        
    return removedNodes


def download_king():
    #king download
    kingDownload = 'http://people.virginia.edu/~wc9c/KING/Linux-king.tar.gz'

    king = shutil.which('king')
    if king is None:
  
        print('King not present, downloading ...')
        cmd = 'wget -O king.tar.gz ' + kingDownload
        subprocess.call(shlex.split(cmd))
        
        print('Extracting...')
        cmd = 'tar -xzvf king.tar.gz -C . '
        subprocess.call(shlex.split(cmd))
        cmd =  'mv king /usr/local/bin'
        
        try:
            subprocess.call(shlex.split(cmd))
        except:
            cmd = 'sudo ' + cmd
            subprocess.call(shlex.split(cmd))
        cmd = 'rm  ./king.tar.gz '
        subprocess.call(shlex.split(cmd))
        print('Done.')
        
        subprocess.call('king')
 
def return_degree(kinship_value):

    for degree in kinship_range:
        degree_extremes = kinship_range[degree]
        if degree_extremes[0]  <= kinship_value <= degree_extremes[1]:
            return degree
    return 4



def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c).copy()
