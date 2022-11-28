from utils import mapcount,subprocess,make_sure_path_exists,write_fam_samplelist,tmp_bash
import networkx as nx
import numpy as np
import os.path,shlex,shutil

def kinship(args):
    '''
    Builds a new plink file for king and performs king analysis
    Returns:
    args.kin (if missings) : output of KING --related
    args.duplicates : list of duplicates
    args.related_couples : tsv file with related couples
    args.related_individuals : list of related individuals
    '''
    
    # Return LIST OF RELATED INDIVIDUALS TO BE REMOVED AND LIST OF DUPLICATES
    args.related_couples = os.path.join(args.kinPath,args.name  +'_related_couples_' + str(args.degree) + '.txt')
    args.duplicates = os.path.join(args.kinPath,args.name  +'_duplicates.txt')
    if not os.path.isfile(args.related_couples) or not os.path.isfile(args.duplicates) or args.force:
        if not args.kin:
            args.kin = os.path.join(args.kinPath,args.name + '.kin0')          
       
        if not os.path.isfile(args.kin):
            print('generating kinship relations file to ',args.kin)
            cmd = f'king -b {args.bed}  --related --degree 2  --cpus {args.cpus} --prefix {args.kin.split(".kin0")[0]}'
            subprocess.call(shlex.split(cmd))
                
        deg_cmd = f"cat  {args.kin}  | grep -vw 3rd |cut -f 2,4| sed -E 1d  > {args.related_couples}"
        tmp_bash(deg_cmd)
        # cuts IID of duplicates, returns uniques and generates duplicate.fam file
        dup_cmd = f"""cat {args.kin} | grep Dup/MZ | cut -f 2,4 |grep -o -E '\\w+' | sort -u -f >""" + args.duplicates
        print(dup_cmd)
        tmp_bash(dup_cmd)
    else:
        args.logging.info("related info already generated")
    
    print(f'Degree {args.degree} related couples : {mapcount(args.related_couples)}')
    print(f"Total duplicates : {mapcount(args.duplicates)}")

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
        print(f'removing non finns before starting to trim nodes {args.all_outliers} {mapcount(args.all_outliers)}')
        remove_nodelist = np.loadtxt(args.all_outliers,dtype = str,usecols = 0)
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
        np.savetxt(args.related_individuals,min(related_greedy,related_nodes),fmt = '%s' )


    else:
        args.logging.info('list of related samples already generated.')
        
    print(f'Degree {args.degree} related individuals : {mapcount(args.related_individuals)}')
    

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
    edges = g.number_of_edges() 
    while edges >0:
        #find highest degree node
        maxNode = max(degrees, key=degrees.get)
        removedNodes.append(maxNode)
        
        for neighbor in g[maxNode]:
            degrees[neighbor] -= 1
            
        #delete node from degree dict and from network
        del degrees[maxNode]
        g.remove_node(maxNode)
        
        # update number of edges
        edges =g.number_of_edges()
        
    return removedNodes

def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c).copy()
