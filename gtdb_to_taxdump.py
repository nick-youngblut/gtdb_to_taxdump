#!/usr/bin/env python
from __future__ import print_function
import sys,os
import argparse
import logging
import csv
import urllib.request
import codecs


desc = 'Converting GTDB taxonomy to NCBI taxdump format'
epi = """DESCRIPTION:
Convert Genome Taxonomy Database (GTDB) taxonomy files
to NCBI taxdump format (names.dmp & nodes.dmp).

The input can be >=1 tsv file with the taxonomy
or >=1 url to the tsv file(s).

The *.dmp files are written to `--outdir`.
A tab-delim table of taxID info is written to STDOUT.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('tax_file', metavar='tax_file', type=str, nargs='+',
                    help='>=1 taxonomy file (or url)')
parser.add_argument('-o', '--outdir', type=str, default='.',
                    help='Output directory (Default: %(default)s)')
parser.add_argument('-e', '--embl-code', type=str, default='XX',
                    help='embl code to use for all nodes (Default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)



class Graph(object):
    def __init__(self, graph_dict=None):
        """ initializes a graph object 
            If no dictionary or None is given, 
            an empty dictionary will be used
        """
        self.__ranks = {'d__' : 'superkingdom',
                        'p__' : 'phylum',
                        'c__' : 'class',
                        'o__' : 'order',
                        'f__' : 'family',
                        'g__' : 'genus',
                        's__' : 'species'}

        if graph_dict == None:
            graph_dict = {}
        self.__graph_dict = graph_dict
        self.__graph_nodeIDs = {}
        self.__seen = {}

    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []
            self.__graph_nodeIDs[vertex] = len(self.__graph_nodeIDs.keys())+1

    def add_edge(self, vertex1, vertex2):
        """ assumes that edge is of type set, tuple or list; 
            between two vertices can be multiple edges! 
        """
        try:
            self.__graph_dict[vertex1].append(vertex2)
        except KeyError:
            self.__graph_dict[vertex1] = [vertex2]
            
    def __generate_edges(self):
        """ A static method generating the edges of the 
            graph "graph". Edges are represented as sets 
            with one (a loop back to the vertex) or two 
            vertices 
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                 if {neighbour, vertex} not in edges:
                    edges.append({vertex, neighbour})
        return edges

    def __str__(self):
        res = "vertices: "
        for k in self.__graph_dict.keys():
            res += str(k) + " "
        res += "\nvertex UIDs: "
        for k in self.__graph_dict:
            res += str(self.__graph_nodeIDs[k]) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res

    def get_rank(self, vertex):
        """ Getting rank based on GTDB prefixes """
        for r in self.__ranks:
            if vertex.startswith(r):
                return self.__ranks[r]
        return 'subspecies'

    def iter_graph(self, vertex):
        """ General iteration of all nodes in the graph """
        if vertex == 'root':
            self.__seen = {}
        for child in self.__graph_dict[vertex]:
            if child not in self.__seen:
                print('Parent: {}; Child: {}'.format(vertex, child))
                self.iter_graph(child)
                self.__seen[child] = 1
    
    def _write_dmp_iter(self, vertex, outName, outNode, embl_code='XX'):
        for child in self.__graph_dict[vertex]:
            if child in self.__seen:
                continue
            self.__seen[child] = 1
            # names
            names = [str(self.__graph_nodeIDs[child]), child, '', 'scientific name']
            outName.write('\t|\t'.join(names) + '\t|\n')
            # nodes
            rank = self.get_rank(child)
            nodes = [self.__graph_nodeIDs[child], self.__graph_nodeIDs[vertex],
                     rank, embl_code, 0, 0, 11, 1, 1, 0, 0, 0]
            outNode.write('\t|\t'.join([str(x) for x in nodes]) + '\t|\n')
            # children
            self._write_dmp_iter(child, outName, outNode, embl_code)
            
    def write_dmp(self, outdir='.', embl_code='XX'):   
        """ Writing names.dmp & nodes.dmp """
        names_file = os.path.join(outdir, 'names.dmp')
        nodes_file = os.path.join(outdir, 'nodes.dmp')
        with open(names_file, 'w') as outName, open(nodes_file, 'w') as outNode:
            # iterating over all vertices starting at the root
            ## writing root
            ### names
            names = [str(self.__graph_nodeIDs['root']), 'all', '', 'synonym']
            outName.write('\t|\t'.join(names) + '\t|\n')
            names = [str(self.__graph_nodeIDs['root']), 'root', '', 'scientific name']
            outName.write('\t|\t'.join(names) + '\t|\n')
            ### nodes
            nodes = [self.__graph_nodeIDs['root'], 1, 'no rank', embl_code,
                     0, 0, 11, 1, 1, 0, 0, 0]
            outNode.write('\t|\t'.join([str(x) for x in nodes]) + '\t|\n')
            ## other nodes
            self._write_dmp_iter('root', outName, outNode, embl_code)
        return names_file, nodes_file

    def _to_tbl_iter(self, vertex):        
        for child in self.__graph_dict[vertex]:
            if child in self.__seen:
                continue
            self.__seen[child] = 1
            # tbl row
            x = [str(self.__graph_nodeIDs[child]), child, self.get_rank(child)]
            print('\t'.join(x))
            # children
            self._to_tbl_iter(child)
    
    def to_tbl(self):
        """ Writing table of values [taxID, name, rank] """
        ## writing header
        x = ['taxID', 'name', 'rank']
        print('\t'.join(x))
        ## writing root
        self.__seen = {}
        x = [str(self.__graph_nodeIDs['root']), 'root', 'no rank']
        print('\t'.join(x))
        self._to_tbl_iter('root')
            
def find_all_paths(self, start_vertex, end_vertex, path=[]):
    """ find all paths from start_vertex to 
        end_vertex in graph """
    graph = self.__graph_dict 
    path = path + [start_vertex]
    if start_vertex == end_vertex:
        return [path]
    if start_vertex not in graph:
        return []
    paths = []
    for vertex in graph[start_vertex]:
        if vertex not in path:
            extended_paths = self.find_all_paths(vertex, 
                                                 end_vertex, 
                                                 path)
            for p in extended_paths: 
                paths.append(p)
    return paths

def load_gtdb_tax(infile, graph):
    """ loading gtdb taxonomy & adding to DAG """

    try:
        ftpstream = urllib.request.urlopen(infile)
        inF = csv.reader(codecs.iterdecode(ftpstream, 'utf-8'))
    except ValueError:
        inF = open(infile)
        
    for i,line in enumerate(inF):
        # parsing
        try:
            line = line.rstrip()
        except AttributeError:
            line = line[0].rstrip()
        if line == '':
            continue
        line = line.split('\t')
        if len(line) < 2:
            msg = 'Line{} does not contain >=2 columns'
            raise ValueError(msg.format(i))            
        tax = line[1].split(';')
        if len(tax) < 7:
            msg = 'WARNING: Line{}: taxonomy length is <7'
            logging.info(msg.format(i))
        tax.append(line[0])
        # adding taxonomy to graph
        for ii,cls in enumerate(tax):
            graph.add_vertex(cls)
            if ii == 0:
                graph.add_edge('root', cls)
            else:
                graph.add_edge(tax[ii-1], cls)
    try:
        inF.close()
    except AttributeError:
        pass
                    
def main(args):
    # creating DAG
    graph = Graph()
    graph.add_vertex('root')
    for F in args.tax_file:
        logging.info('Loading: {}'.format(F))
        load_gtdb_tax(F, graph)
    # transversing DAG and writing node.dmp & names.dmp
    names_file, nodes_file = graph.write_dmp(args.outdir,
                                             embl_code=args.embl_code)
    logging.info('File written: {}'.format(names_file))
    logging.info('File written: {}'.format(nodes_file))
    # writing standard table of taxIDs
    graph.to_tbl()
    
         
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)


# node.dmp ("\t|\t" delimited)
#tax_id => nodeID              
#parent tax_id => parent_nodeID
#rank => steps from root (no rank, domain, phylum, class, order, family, genus, species, strain)
#embl code => XX?			
#division id => 0
#inherited div flag => 0
#genetic code id => 11
#inherited GC  flag => 1 (except for root)
#mitochondrial genetic code id => 1
#inherited MGC flag  => 0
#GenBank hidden flag => 0
#hidden subtree root flag => 0
#comments => ''

# names.dmp
#tax_id => nodeID
#name_txt => taxonomic classification
#unique name => taxonomic classification (full lineage)
#name class => "scientific name"
## Notes:
# * need "all"
