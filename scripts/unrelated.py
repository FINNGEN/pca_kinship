import networkx as nx
import os,argparse
from utils import progressBar,file_exists


def filter_kinship(out_dir,kin_file, kinship_filter):
    with open(kin_file) as i:header=i.readline().strip().split()
    idx = [header.index(elem) for elem in ['ID1','ID2','Kinship']]
    print(list(zip(idx,[header[i] for i in idx])))
    
    out_file = os.path.join(out_dir,f"unrelated_pairs_{kinship_filter}.txt")
    if not  os.path.isfile(out_file):
        with open(out_file,'wt') as o,open(kin_file) as i:
            next(i)
            for line in i:
                data = line.strip().split()
                s1,s2,kinship= [data[i] for i in idx]
                if float(kinship) > kinship_filter:
                    o.write('\t'.join((s1,s2,kinship))+'\n')
    print(out_file)
    return out_file


def return_related(edgelist):

    g = nx.read_weighted_edgelist(edgelist)
    print(f"edges:{g.number_of_edges()}")
    # native nx algorithm vs greedy algorithm
    unrelated_nodes = []
    related_greedy = []
    components = list(connected_component_subgraphs(g))
    n_components = len(components)
    print(f"components:{n_components}")
    for i,subgraph in enumerate(components):
        progressBar(i+1, n_components, bar_length=20)
        unrelated_nodes += nx.maximal_independent_set(subgraph)
        related_greedy += greedy_algorithm(subgraph)
    print('\nDone.')
    #sanity checks
    sanity_check(g,unrelated_nodes)
    sanity_check(g,g.nodes() - related_greedy)
    related_nodes = list(g.nodes() - set(unrelated_nodes))
    print(f'{len(related_nodes)} native related')
    print(f'{len(related_greedy)} greedy related')

    final_related = min([related_nodes,related_greedy], key=len)
    print(f'{len(final_related)} final related')
    return final_related
    
def sanity_check(graph,nodes):
    assert graph.subgraph(nodes).number_of_edges() == 0

    
def greedy_algorithm(g):
    """
    Removes sequentially node with highest degree until there are no nodes left
    """
    
    degrees = dict(g.degree())
    removedNodes = []
    #edges = 
    while g.number_of_edges()   >0:
        #find highest degree node
        maxNode = max(degrees, key=degrees.get)
        removedNodes.append(maxNode)
        for neighbor in g[maxNode]:degrees[neighbor] -= 1         
        #delete node from degree dict and from network
        del degrees[maxNode]
        g.remove_node(maxNode)
    return removedNodes

def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c).copy()

def main(args):

    edgelist = filter_kinship(args.out_dir,args.kinship,args.kin_filter)
    related = return_related(edgelist)
    out_file = os.path.join(args.out_dir,f"related_samples_{args.kin_filter}.txt")
    with open(out_file,'wt') as f:
        for sample in related:f.write(sample + '\n')

if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description="unrelated pipeline.")
    parser.add_argument("--kinship", type=file_exists, help =  "Path to csv file with sample,batch", default ='/home/pete/r12/kinship/release/data/finngen_R12.kin0' )
    parser.add_argument('--kin-filter',type=float,help='Degree for Kinship',default = 0.25)
    parser.add_argument('-o',"--out_dir",type = str, help = "Folder in which to save the results", required = True)

    args = parser.parse_args()
    main(args)
