#!/usr/bin/env python
from __future__ import print_function
# batteries
import os
import sys
import re
import shutil
import tarfile
import argparse
import logging
import random
import csv
import urllib.request
import codecs
import functools
import multiprocessing as mp
from collections import Counter
# 3rd party
import networkx as nx
from networkx.algorithms.dag import descendants
from networkx.algorithms.lowest_common_ancestors import lowest_common_ancestor
from networkx.algorithms.shortest_paths.unweighted import bidirectional_shortest_path

desc = 'Mapping between NCBI & GTDB taxonomies'
epi = """DESCRIPTION:
Using the GTDB metadata table (which contains both NCBI and GTDB taxonomies)
to map taxonomic classifications between the 2 taxonomies.

For example, determine GTDB equivalent of the NCBI taxonomic classifications:
* Bacillus
* Xanthomonas oryzae
* Burkholderiaceae

Algorithm:
* Read in GTDB metadata (bac and/or arc)
* Convert NCBI & GTDB taxonomies of each entry to a directed acyclic graph (DAG)
  * "accession" used for iips of each taxonomy tree
* Read in taxonomy queries (a file w/ 1 entry per line)
* For each query:
  * Found in reference taxonomy (hit/no-hit)?
  * If in reference taxonomy, get all tips (all accessions) for that node in the DAG
  * For all tips, get LCA in target taxonomy (if NCBI queries, then GTDB is target)
    * Only <= (--max-tips) used to save on time. Tips are subsampled randomly.
    * "fuzzy" LCAs allowed, in which only >= (--fraction) of the tips must have that LCA
  * If no match in reference taxonomy or no LCA that is >= (--fraction), then the 
    target LCA is "unclassified"

Output table columns:
  * ncbi_taxonomy
    * NCBI taxonomy name
  * gtdb_taxonomy
    * GTDB taxonomy name
  * lca_frac
    * Fraction of tips with that LCA
  * target_tax_level
    * The taxonomic level of the target (eg., genus or species)
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('tax_queries', metavar='tax_queries', type=str,
                    help='List of taxa to query (1 per line)')
parser.add_argument('gtdb_metadata', metavar='gtdb_metadata', type=str, nargs='+',
                    help='>=1 gtdb-metadata file (or url)')
parser.add_argument('-q', '--query-taxonomy', type=str, default='ncbi_taxonomy',
                    choices=['ncbi_taxonomy', 'gtdb_taxonomy'],
                    help='Taxonomy of the query list (Default: %(default)s)')
parser.add_argument('-f', '--fraction', type=float, default=0.9,
                    help='Homogeneity of LCA (fraction) in order to be used (Default: %(default)s)')
parser.add_argument('-m', '--max-tips', type=int, default=100,
                    help='Max no. of tips used for LCA determination. If more, subsampling w/out replacement (Default: %(default)s)')
parser.add_argument('-p', '--procs', type=int, default=1,
                    help='No. of parallel processes (Default: %(default)s)')
parser.add_argument('-v', '--verbose', action='store_true', default=False,
                    help='Verbose output (Default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def format_taxonomy(T, hierarchy, acc):
    """
    Formatting taxonomy to conform to a set hierarchy
    """
    Tx = ['' for i in range(len(hierarchy))]
    for i,x in enumerate(hierarchy[:-1]):
        if len(T) < i + 1 or T[i] == '' or T[i] == 'unclassified':
            Tx[i] = ':'.join([x, acc])
        else:
            Tx[i] = T[i]
    Tx[-1] = acc
    return Tx

def add_taxonomy(line, line_num, header, G, tax='ncbi_taxonomy'):
    """
    Adding taxonomy nodes/edits to the graph
    """
    regex = re.compile(r'^[dpcofgs]__')
    hierarchy = ['domain', 'phylum', 'class', 'order',
                 'family', 'genus', 'species', 'strain']
    # checking taxonomy format
    acc = line[header['accession']]
    T = line[header[tax]].split(';')
    T = [regex.sub('', x) for x in T]
    T = format_taxonomy(T, hierarchy, acc)
    # adding taxonomy to graph
    for i in range(len(hierarchy)):
        # adding node
        G[tax].add_node(T[i], taxonomy=hierarchy[i])
        # adding edge
        if i == 0:
            G[tax].add_edge('root', T[i])
        else:
            G[tax].add_edge(T[i-1], T[i])

def dl_uncomp(url):
    """
    Downloading and extracting GTDB metadata tarball
    """
    file_tmp = urllib.request.urlretrieve(url, filename=None)[0]
    tmpdir = os.path.dirname(file_tmp)
    outdir = os.path.basename(url)
    outdir = re.sub(r'\.tar.gz$', '', outdir)
    outdir = os.path.join(tmpdir, outdir)
    tar = tarfile.open(file_tmp)
    tar.extract(tar.getmembers()[0], outdir)
    inF = open(os.path.join(outdir, tar.getnames()[0]))
    return inF,outdir
            
def load_gtdb_metadata(infile, G):
    """
    Loading gtdb taxonomy & adding to DAG 
    """
    # input as file or url
    try:
        inF,tmpdir = dl_uncomp(infile)
    except ValueError:
        inF = open(infile)
        tmpdir = None
    # reading
    header = {}
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
        # header
        if i == 0:
            header = {x:ii for ii,x in enumerate(line)}
            continue
        # Adding taxonomies to graphs
        add_taxonomy(line, i, header, G, tax='gtdb_taxonomy')
        add_taxonomy(line, i, header, G, tax='ncbi_taxonomy')
    # closing
    try:
        inF.close()
    except AttributeError:
        pass
    if tmpdir is not None and os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    return G

def DiGraph_w_root():
    """
    directed graph with a root node
    """
    G = nx.DiGraph()
    G.add_node('root')
    return G

def lca_frac_pass(D, lca_frac):
    """
    Determine which, if any, of the LCAs pass the homogeneity cutoff.
    "homogeneity" means the fraction of tips with the target LCA 
      (eg., 90% of all tips have this LCA).
    If the cutoff is not passed, then returning [None,None]
    """
    D = Counter(D)
    try:
        mc = D.most_common(1)
    except IndexError:
        return [None,None]
    try:
        frac = mc[0][1] / float(sum(D.values()))
    except IndexError:
        return [None,None]
    if frac >= lca_frac:
        return [mc[0][0], str(round(frac, 3))]
    else:
        return [None,None]

def lca_many_nodes(G, nodes, lca_frac=1.0):
    """
    algorithm: using closest distance to 'root'
    """
    hierarchy = ['root', 'domain', 'phylum', 'class', 'order',
                 'family', 'genus', 'species', 'strain']    
    T = [{} for x in hierarchy]
    # getting path from nodes to root
    for n in nodes:
        path = bidirectional_shortest_path(G, 'root', n)
        for i,node in enumerate(path):
            try:
                T[i][node] += 1
            except KeyError:
                T[i][node] = 1
                
    # from tip to root, which passess LCA cutoff?
    ## note: species is lowest possible level
    for i in range(len(T)-1)[::-1]:
        lca = lca_frac_pass(T[i], lca_frac)
        if lca[0] is not None:
            return lca + [hierarchy[i]]
    raise ValueError('Cannot find LCA for nodes: {}'.format(','.join(nodes)))

def _query_tax(tax_queries, G, qtax, ttax, lca_frac=1.0, max_tips=100, verbose=False):
    """
    Querying list of taxonomic names.
    """
    pid = os.getpid()
    idx = {}
    status = {'hit' : 0, 'no hit': 0}
    # iterating queries
    for Q in tax_queries:
        tips = []
        try:
            # getting descendents of the node
            tips = [desc for desc in descendants(G[qtax], Q) if \
                    G[qtax].nodes[desc]['taxonomy'] == 'strain']
            status['hit'] += 1
        except nx.exception.NetworkXError:
            status['no hit'] += 1
        # if tips, getting LCA in target-taxonomy
        n_tips = len(tips)
        if n_tips > 0:
            if n_tips > max_tips:                    
                tips = random.sample(tips, k=max_tips)
            LCA = lca_many_nodes(G[ttax], tips, lca_frac=lca_frac)
            idx[Q] = LCA
        else:
            idx[Q] = ['unclassified', 'NA', 'NA']
        # status
        x = status['hit'] + status['no hit']
        if verbose and x % 1000 == 0:
            frac = round(float(x) / len(tax_queries) * 100, 2)
            logging.info('PID{}: Queries processed: {} ({}%)'.format(pid, x, frac))       
    # status
    msg = 'PID{}: Finished! Queries={}, Hits={}, No-Hits={}'
    logging.info(msg.format(pid, status['hit'] + status['no hit'],
                            status['hit'], status['no hit']))
    # return
    return idx

def query_tax(tax_queries, G, tax, lca_frac=1.0, max_tips=100, procs=1, verbose=False):
    ttax = 'ncbi_taxonomy' if tax == 'gtdb_taxonomy' else 'gtdb_taxonomy'
    # loading & batching queries
    logging.info('Reading in queries: {}'.format(tax_queries))    
    q_batch = [[] for i in range(procs)]
    with open(tax_queries) as inF:
        for i,line in enumerate(inF):
            line = line.rstrip()
            if line == '' or line == 'root':
                continue
            q_batch[i % procs].append(line)
            # debug
            #if i > 10000:
            #    break
    # query graphs
    logging.info('Querying taxonomies...')
    func = functools.partial(_query_tax, G=G, qtax=tax, ttax=ttax,
                             lca_frac = lca_frac, max_tips=max_tips,
                             verbose=verbose)
    if procs > 1:
        P = mp.Pool(procs)
        idx = P.map(func, q_batch)
    else:
        idx = map(func, q_batch)    
    return idx

def write_table(idx, qtax):
    """
    Writing tab-delim table of taxonomy mappings to STDOUT
    """
    logging.info('Writing table to STDOUT...')
    ttax = 'ncbi_taxonomy' if qtax == 'gtdb_taxonomy' else 'gtdb_taxonomy'    
    print('\t'.join([qtax, ttax, 'lca_frac', 'target_tax_level']))
    for x in idx:
        for k,v in x.items():
            print('\t'.join([k] + v))

def main(args):
    """
    Main interface
    """
    # loading the graphs
    G = {'ncbi_taxonomy' : DiGraph_w_root(),
         'gtdb_taxonomy' : DiGraph_w_root()}
    for F in args.gtdb_metadata:
        logging.info('Loading: {}'.format(F))
        load_gtdb_metadata(F, G)
    # querying
    idx = query_tax(args.tax_queries, G,
                    tax=args.query_taxonomy,
                    lca_frac = args.fraction,
                    max_tips = args.max_tips,
                    procs = args.procs,
                    verbose = args.verbose)
    # writing results
    write_table(idx, qtax=args.query_taxonomy)
             
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

