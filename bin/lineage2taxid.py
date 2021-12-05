#!/usr/bin/env python
# import
## batteries
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging
## 3rd party
import networkx as nx
from networkx.algorithms.dag import descendants
from networkx.algorithms.lowest_common_ancestors import lowest_common_ancestor
from networkx.algorithms.shortest_paths.unweighted import bidirectional_shortest_path
## package
import gtdb_to_taxdump as gtdb2td

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

# argparse
desc = 'Map taxonomic lineages to taxids'
epi = """DESCRIPTION:
If you have a set of taxonomic lineages (eg., generated from GTDB-Tk)
and need then taxids for the taxonomic lineages, then this script is for you!

Example lineages: 
* d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Ruminococcus_E;s__Ruminococcus_E bromii_B
* d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola plebeius_A

Method:
* The nodes.dmp & names.dmp files are loaded as a graph
* For each lineage in the "table_file":
  * From the finest-to-coarsest taxonomic resolution:
    * Find node in graph with that taxonomic classification 
      * Note: captilization invariant
      * If only 1 node (so no duplicate naming), then return that taxid & rank
      * If not, go to next-finest taxonomic resolution
    * If no taxonomic classification can be mapped:
      * Use "NA" and warn the user

Notes:
* The input table file can be compressed via gzip or bzip2
* The input table is assumed to be tab-delimited
* The output is writted to STDOUT
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('table_file', metavar='table_file', type=str,
                    help='Tab-delim input table containing GTDB taxonomic lineages')
parser.add_argument('names_dmp', metavar='names_dmp',
                    help='NCBI names.dmp file. Only needed if providing NCBI taxids')
parser.add_argument('nodes_dmp', type=str, default=None,
                    help='NCBI nodes.dmp file. Only needed if providing NCBI taxids')
parser.add_argument('--lineage-column', type=str, default='classification',
                    help='Column name that contains the lineages')
parser.add_argument('--taxid-column', type=str, default='taxid',
                    help='Name of taxid column that will be appended to the input table')
parser.add_argument('--taxid-rank-column', type=str, default='taxid_rank',
                    help='Name of taxid-rank column that will be appended to the input table')
parser.add_argument('--version', action='version', version='0.0.1')

# functions
def load_dmp(names_dmp_file, nodes_dmp_file):
    """
    Loading NCBI names/nodes dmp files as DAG
    Arguments:
      names_dmp_file : str, names.dmp file
      nodes_dmp_file : str, nodes.dmp file 
    Return:
      network.DiGraph object
    """
    regex = re.compile(r'\t\|\t')
    # nodes
    logging.info('Loading file: {}'.format(names_dmp_file))
    idx = {}    # {taxid : name}
    with gtdb2td.Utils.Open(names_dmp_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            line = regex.split(line)
            idx[int(line[0])] = line[1].lower()
    # names
    logging.info('Loading file: {}'.format(nodes_dmp_file))
    G = nx.DiGraph()
    G.add_node(0, rank = 'root', name = 'root')
    with gtdb2td.Utils.Open(nodes_dmp_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            line = regex.split(line)
            taxid_child = int(line[0])
            taxid_parent = int(line[1])
            rank_child = line[2]
            name_child = idx[taxid_child].lower()
            name_parent = idx[taxid_parent].lower()
            # adding node
            G.add_node(taxid_child, rank=rank_child, name=name_child)
            # adding edge
            if taxid_parent == 1:
                G.add_edge(0, taxid_child)
            else:
                G.add_edge(taxid_parent, taxid_child)
    idx.clear()
    logging.info('  No. of nodes: {}'.format(G.number_of_nodes()))
    logging.info('  No. of edges: {}'.format(G.number_of_edges()))
    return G

def lineage2taxid(lineage, G):    
    lineage = lineage.split(';')    
    for cls in lineage[::-1]:
        cls = cls.lower()
        nodes = [x for x,y in G.nodes(data=True) if y['name'] == cls]
        if len(nodes) == 1:
            return [nodes[0], G.nodes[nodes[0]]['rank']]
    msg = 'Could not find a taxid for lineage: {}'
    logging.warning(msg.format(';'.join(lineage)))
    return ['NA', 'NA']            

def parse_lineage_table(table_file, lineage_column, G,
                        taxid_column, taxid_rank_column):
    """
    Parsing lineage and finding taxid
    """
    logging.info('Parsing file: {}'.format(table_file))
    header = {}
    with gtdb2td.Utils.Open(table_file) as inF:
        for i,line in enumerate(inF):
            line = gtdb2td.Utils.Decode(line).rstrip().split('\t')
            # header
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                try:
                    _ = header[lineage_column]
                except KeyError:
                    msg = 'Cannot find column: {}'
                    raise KeyError(msg.format(lineage_column))
                print('\t'.join(line + [taxid_column, taxid_rank_column]))
                continue
            # body
            lineage = line[header[lineage_column]]
            taxid,rank = lineage2taxid(lineage, G)
            print('\t'.join(line + [str(taxid), str(rank)]))
            # status
            if i > 0 and (i+1) % 100 == 0:
                logging.info('  Records processed: {}'.format(i+1))
            
def main(args):
    """
    Main interface
    """
    # loading dmp as DAG
    G = load_dmp(args.names_dmp, args.nodes_dmp)
    # lineage2taxid
    parse_lineage_table(args.table_file, args.lineage_column, G=G,
                        taxid_column = args.taxid_column,
                        taxid_rank_column = args.taxid_rank_column)
        
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
