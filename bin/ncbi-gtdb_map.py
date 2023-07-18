#!/usr/bin/env python
from __future__ import print_function
# import
## batteries
import os
import re
import shutil
import tarfile
import argparse
import logging
import random
import urllib.request
import functools
import multiprocessing as mp
from collections import Counter
## 3rd party
import networkx as nx
from networkx.algorithms.dag import descendants
from networkx.algorithms.dag import ancestors
from networkx.algorithms.lowest_common_ancestors import lowest_common_ancestor
from networkx.algorithms.shortest_paths.unweighted import bidirectional_shortest_path
from bin import __version__
## package
import gtdb2td

# argparse
desc = 'Map between NCBI & GTDB taxonomies'
epi = """DESCRIPTION:
Using the GTDB metadata table (which contains both NCBI and GTDB taxonomies)
to map taxonomic classifications between the 2 taxonomies.

For example, determine GTDB equivalent of the NCBI taxonomic classifications:
* g__Bacillus
* s__Xanthomonas oryzae
* f__Burkholderiaceae
* s__Methanobrevibacter smithii

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

Notes:
* The input table of queries should be tab-delimited (select the column with --column).
* Query-target matching is caps-invariant (all converted to lower case for matching)!
* The table can have a header (see --header) and can be compressed via gzip or bzip2.
* See https://github.com/nick-youngblut/gtdb_to_taxdump/tree/master/tests/data/ncbi-gtdb for example queries.

Output:
* Output written to --outdir
* `taxonomy_map_summary.tsv`
  * ncbi_taxonomy
    * NCBI taxonomy name
  * gtdb_taxonomy
    * GTDB taxonomy name
  * lca_frac
    * Fraction of tips with that LCA
  * target_tax_level
    * The taxonomic level of the target (eg., genus or species)
* `queries_renamed.tsv` 
  * if --renamed
  * Note: renaming won't work if using NCBI taxids
  * the format of the output table will match the query table

Examples:

# NCBI => GTDB
python ncbi-gtdb_map.py tests/data/ncbi-gtdb/ncbi_tax_queries.txt https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/ar122_metadata_r95.tar.gz https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz

# GTDB => NCBI
python ncbi-gtdb_map.py -q gtdb_taxonomy tests/data/ncbi-gtdb/gtdb_tax_queries.txt https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/ar122_metadata_r95.tar.gz https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz

# NCBI => GTDB, using NCBI taxids
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -pzxvf taxdump.tar.gz
python ncbi-gtdb_map.py --names-dmp taxdump/names.dmp --nodes-dmp taxdump/nodes.dmp tests/data/ncbi-gtdb/ncbi_taxid_queries.txt https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/ar122_metadata_r95.tar.gz https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz

# NCBI => GTDB; no lineage prefix ([dpcofgs]__)
python ncbi-gtdb_map.py -N tests/data/ncbi-gtdb/ncbi_tax_queries_noPrefix.txt https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/ar122_metadata_r95.tar.gz https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('tax_queries', metavar='tax_queries', type=str,
                    help='List of taxa to query (1 per line)')
parser.add_argument('gtdb_metadata', metavar='gtdb_metadata', type=str, nargs='+',
                    help='>=1 gtdb-metadata file (or url)')
parser.add_argument('-o', '--outdir', type=str, default='ncbi-gtdb',
                    help='Output file directory (Default: %(default)s)')
parser.add_argument('-q', '--query-taxonomy', type=str, default='ncbi_taxonomy',
                    choices=['ncbi_taxonomy', 'gtdb_taxonomy'],
                    help='Taxonomy of the query list (Default: %(default)s)')
parser.add_argument('-f', '--fraction', type=float, default=0.90,
                    help='Homogeneity of LCA (fraction) in order to be used (Default: %(default)s)')
parser.add_argument('-m', '--max-tips', type=int, default=100,
                    help='Max no. of tips used for LCA determination. If more, subsampling w/out replacement (Default: %(default)s)')
parser.add_argument('-c', '--column', type=int, default=1,
                    help='Column containing the queries;' + \
                         ' assuming a tab-delim table (Default: %(default)s)')
parser.add_argument('-H', '--header', action='store_true', default=False,
                    help='Header in table of queries (Default: %(default)s)?')
parser.add_argument('-P', '--prefix', type=str, default='',
                    help='Add prefix to all queries such as "s__" (Default: %(default)s)')
parser.add_argument('-N', '--no-prefix', action='store_true', default=False,
                    help='Strip off any [dpcofgs]__ taxonomic prefixes? (Default: %(default)s)')
parser.add_argument('--completeness', type=float, default=50.0,
                    help='Only include GTDB genomes w/ >=X CheckM completeness (Default: %(default)s)')
parser.add_argument('--contamination', type=float, default=5.0,
                    help='Only include GTDB genomes w/ <X CheckM contamination (Default: %(default)s)')
parser.add_argument('--names-dmp', type=str, default=None,
                    help='NCBI names.dmp file. Only needed if providing NCBI taxids (Default: %(default)s)')
parser.add_argument('--nodes-dmp', type=str, default=None,
                    help='NCBI nodes.dmp file. Only needed if providing NCBI taxids (Default: %(default)s)')
parser.add_argument('-r', '--rename', action='store_true', default=False,
                    help='Write query file with queries renamed? (Default: %(default)s)')
parser.add_argument('-p', '--procs', type=int, default=1,
                    help='No. of parallel processes (Default: %(default)s)')
parser.add_argument('-v', '--verbose', action='store_true', default=False,
                    help='Verbose output (Default: %(default)s)')
parser.add_argument('--version', action='version', version=__version__)
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


# functions

def format_taxonomy(T, hierarchy: list, acc: str, no_prefix=False) -> list:
    """
    Formatting taxonomy to conform to a set hierarchy
    Params:
      T: iterable, taxonomy
      hierarchy: taxonomic hierarchy 
      acc: accession
      no_prefix: strip taxonomic prefixes?
    Return:
      List of taxonomic names in the hierarchy; [taxonomy_level1, taxonomy_level2, ...]
    """
    regex = re.compile(r'^[dpcofgs]__$')
    regex2 = re.compile(r'^X*[dpcofgs]__')    
    Tx = ['' for i in range(len(hierarchy))]
    for i,x in enumerate(hierarchy[:-1]):
        if len(T) < i + 1 or T[i] == '' or T[i] == 'unclassified' \
           or regex.search(T[i]):
            Tx[i] = '__'.join(['X' + x[0], acc])
        else:
            Tx[i] = T[i]
        if no_prefix is True:
            Tx[i] = regex2.sub('', Tx[i])
    Tx[-1] = acc
    return Tx

def add_taxonomy(line: list, line_num: int, header: dict, G, 
                 tax='ncbi_taxonomy', no_prefix=False) -> None:
    """
    Add taxonomy nodes/edits to the graph.
    Params:
      line: iterable, line of metadata file
      line_num: int, line number of metadata file
      header : dict, index of header for metadata file
      G : taxonomy graph
      tax : ncbi or gtdb?
    Return:
      None; in-place edit of G
    """
    hierarchy = ['domain', 'phylum', 'class', 'order',
                 'family', 'genus', 'species', 'strain']
    # check taxonomy format
    acc = line[header['accession']]
    T = line[header[tax]].split(';')
    T = format_taxonomy(T, hierarchy, acc, no_prefix=no_prefix)
    # add taxonomy to graph
    for i in range(len(hierarchy)):        
        # add node
        G[tax].add_node(T[i].lower(), taxonomy=hierarchy[i], orig_name=T[i])
        # add edge
        if i == 0:
            G[tax].add_edge('root', T[i].lower())
        else:
            G[tax].add_edge(T[i-1].lower(), T[i].lower())

def dl_uncomp(url: str) -> tuple:
    """
    Download and extracting GTDB metadata tarball.
    Params:
      url: str, url of GTDB metadata tarball
    Return: 
      Tuple of (metadata_filehandle, output_directory)
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
            
def load_gtdb_metadata(infile: str, G, completeness: float, 
                       contamination: float, no_prefix: bool=False):
    """
    Load gtdb taxonomy & adding to DAG.
    Params:
      infile: input file url or path
      G: graph object
      completeness: min completeness
      contamination: max contamination
      no_prefix: strip taxonomic prefixes?
    Return:
      Graph object (DAG)
    """
    logging.info(f'Loading: {infile}')
    # input as file or url
    try:
        inF,tmpdir = dl_uncomp(infile)
    except ValueError:
        inF = gtdb2td.Utils.Open(infile)
        tmpdir = None
    # reading
    stats = {'passed' : 0, 'completeness' : 0,
             'contamination' : 0, 'no ncbi tax' : 0}
    header = {}
    for i,line in enumerate(inF):        
        # parsing
        line = gtdb2td.Utils.Decode(line)
        try:
            line = line.rstrip()
        except AttributeError:
            line = line[0].rstrip()
        if line == '':
            continue
        line = line.split('\t')
        if len(line) < 2:
            msg = 'Line{} does not contain >=2 columns'
            raise ValueError(msg.format(i+1))
        # header
        if i == 0:
            header = {x:ii for ii,x in enumerate(line)}
            continue
        # filtering out records lacking an NCBI taxonomy
        try:
            X = line[header['ncbi_taxonomy']]
        except KeyError:
            raise KeyError('Cannot find "ncbi_taxonomy"')
        if X == 'none':
            stats['no ncbi tax'] += 1
            continue
        # filtering by checkM stats
        try:
            X = line[header['checkm_completeness']]
        except KeyError:
            raise KeyError('Cannot find "checkm_completeness"')
        if float(X) < completeness:
            stats['completeness'] += 1
            continue
        try:
            X = line[header['checkm_contamination']]
        except KeyError:
            raise KeyError('Cannot find "checkm_contamination"')
        if float(X) >= contamination:
            stats['contamination'] += 1
            continue
        # Adding taxonomies to graphs        
        add_taxonomy(line, i, header, G, tax='gtdb_taxonomy',
                     no_prefix=no_prefix)
        add_taxonomy(line, i, header, G, tax='ncbi_taxonomy',
                     no_prefix=no_prefix)
        stats['passed'] += 1
    # closing
    try:
        inF.close()
    except AttributeError:
        pass
    if tmpdir is not None and os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    # stats
    msg = '  Entries lacking an NCBI taxonomy: {}'
    logging.info(msg.format(stats['no ncbi tax']))
    msg = '  Completeness-filtered entries: {}'
    logging.info(msg.format(stats['completeness']))
    msg = '  Contamination-filtered entries: {}'
    logging.info(msg.format(stats['contamination']))
    logging.info(f'  Entries used: {stats["passed"]}')
    return G

def DiGraph_w_root():
    """
    Create directed  graph with a root node
    """
    G = nx.DiGraph()
    G.add_node('root')
    return G

def lca_frac_pass(D: dict, lca_frac: float) -> None:
    """
    Determine which, if any, of the LCAs pass the homogeneity cutoff.
    "homogeneity" means the fraction of tips with the target LCA 
      (eg., 90% of all tips have this LCA).
    If the cutoff is not passed, then returning [None,None]
    Params:
      D: LCA dict
      lca_frac: fraction of decendents that must have the same taxonomy
    Return:
      None
    """
    D = Counter(D)
    try:
        mc = D.most_common(1)
    except IndexError:
        return [None,None]
    if re.search(r'^[Xx][pcofgs]__', mc[0][0]):
        return [None,None]
    try:
        frac = mc[0][1] / float(sum(D.values()))
    except IndexError:
        return [None,None]
    if frac >= lca_frac:
        return [mc[0][0], str(round(frac, 3))]
    else:
        return [None,None]

def lca_many_nodes(G, nodes: list, lca_frac: float=1.0) -> list:
    """
    Algorithm: using closest distance to 'root'
    Params:
      G: graph object
      nodes: list of nodes
      lca_frac: fraction of decendents that must have the same taxonomy
    Return:
      List of LCAs that pass the homogeneity cutoff
    """
    hierarchy = ['root', 'domain', 'phylum', 'class', 'order',
                 'family', 'genus', 'species', 'strain']    
    T = [{} for x in hierarchy]
    # Get path from nodes to root
    for n in nodes:
        path = bidirectional_shortest_path(G, 'root', n)
        for i,node in enumerate(path):
            try:
                T[i][node] += 1
            except KeyError:
                T[i][node] = 1
    ## From finest to coarsest, does any classification pass the lca_frac?
    ### note: species is lowest possible level
    for i in range(len(T)-1)[::-1]:
        lca = lca_frac_pass(T[i], lca_frac)
        if lca[0] is not None and lca[0] != 'root':
            try:
                lca[0] = G.nodes[lca[0]]['orig_name']
            except KeyError:
                msg = 'Cannot find "orig_name" for "{}"'
                raise KeyError(msg.format(lca[0]))
            lineage = [[y for y in x.keys()][0] for x in T[:i+1]]
            return lca + [hierarchy[i], ';'.join(lineage[1:])]
    logging.warning('Cannot find LCA for nodes: {}'.format(','.join(nodes)))
    return ['unclassified', 'NA', 'NA', 'NA']

def _query_tax(tax_queries: list, G, qtax: str, ttax: str, 
               lca_frac: float=1.0, max_tips: int=100,
               verbose: bool=False) -> dict:
    """
    Querying list of taxonomic names.
    Params:
      tax_queries: list of taxonomic names to query
      G: graph containing taxonomy
      qtax: query taxonomy column name
      ttax: target taxonomy column name
      lca_frac: fraction of decendents that must have the same taxonomy
      max_tips: maximum number of tips to consider
      verbose: print progress updates?
    Return:
      {query : [new_taxname, lca_score, rank]}
    """
    pid = os.getpid()
    idx = {}
    status = {'hit' : 0, 'no hit': 0}
    # iterating queries
    for Q in tax_queries:
        tips = []
        try:            
            # getting descendents of the node
            tips = [desc for desc in descendants(G[qtax], Q[0]) if \
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
            idx[Q[1]] = LCA
        else:
            idx[Q[1]] = ['unclassified', 'NA', 'NA', 'NA']
        # status
        x = status['hit'] + status['no hit']
        if verbose and x % 1000 == 0:
            frac = round(float(x) / len(tax_queries) * 100, 2)
            msg = 'PID{}: Queries processed: {} ({}%)'
            logging.info(msg.format(pid, x, frac))       
    # status
    msg = 'PID{}: Finished! Queries={}, Hits={}, No-Hits={}'
    logging.info(msg.format(pid, status['hit'] + status['no hit'],
                            status['hit'], status['no hit']))
    # return
    return idx

def queries_taxid2species(queries: dict, tax_graph) -> dict:
    """
    If queries were taxids, getting taxonomic classification for each query.
    Params:
      queries: dict of {query_id: taxid}
      tax_graph: networkx graph of taxonomy
    Return:
      {query_id: [taxname, rank]}
    """
    logging.info('Converting query taxids to species-level classifications')
    queries_new = {}
    for q,cnt in queries.items():
        try:
            q = int(list(q)[0])
        except ValueError:
            msg = 'Cannot convert "{}" to an integer. Is it a taxid?'
            raise ValueError(msg.format(list(q)[0]))
        try:
            node = tax_graph.nodes[q]
        except KeyError:
            msg = 'Cannot find "{}" in NCBI taxdump graph'
            logging.warning(msg.format(q))
            continue
        # if species level, getting species level classification
        try:
            rank = node['rank']
        except KeyError:
            rank = ''
            logging.warning('Cannot find rank for {}'.format(q))            
        if rank == 'species':
            queries_new[(node['name'].lower(), node['name'])] = cnt        
    queries.clear()
    logging.info('  No. of queries: {}'.format(sum(queries_new.values())))
    logging.info('  No. of de-rep queries: {}'.format(len(queries_new.keys())))
    return queries_new
                
def query_tax(tax_queries: str, G, tax: str, lca_frac: float=1.0, 
              max_tips: int=100, column: int=1, header: bool=False, 
              tax_graph=None, prefix: str='',
              procs: int=1, verbose: bool=False) -> dict:
    """
    Querying list of taxonomic names.   
    Params:
      tax_queries: text file of taxonomic queries
      G: nx.graph, graphs for each of the 2 taxonomies
      tax: either 'ncbi_taxonomy' or 'gtdb_taxonomy'
      lca_frac: taxonomic "purity" of node required to consider it the LCA
      max_tips: the max number of tips considered for determining the LCA
      column: the column in the tax_queries file that contains the taxonomic classifications
      header: does the tax_queries file contain a header?
      tax_graph: a graph that will be used for converting NCBI taxids to species-level queries
      prefix: add a prefix to all of the queries (eg., "s__")?      
    Return:
      {(node_name_lowercase, node_name) : count}   
    """
    ttax = 'ncbi_taxonomy' if tax == 'gtdb_taxonomy' else 'gtdb_taxonomy'
    # loading & batching queries
    logging.info(f'Reading in queries: {tax_queries}')
    n_queries = 0
    queries = {}
    with gtdb2td.Utils.Open(tax_queries) as inF:
        for i,line in enumerate(inF):
            if i == 0 and header:
                continue
            line = gtdb2td.Utils.Decode(line)
            line = line.rstrip().split('\t')[column - 1]
            if line == '' or line == 'root':
                continue
            line = (prefix + line.lower(), line)
            try:
                queries[line] += 1
            except KeyError:
                queries[line] = 1
            # debug
            #if i > 10000:
            #    break
    logging.info('No. of queries: {}'.format(sum(queries.values())))
    logging.info('No. of de-rep queries: {}'.format(len(queries.keys())))
    # converting to species-level queries if tax_graph
    if tax_graph is not None:
        queries = queries_taxid2species(queries, tax_graph)
    # batching
    logging.info('Batching queries...')
    q_batch = [[] for i in range(procs)]
    for i,q in enumerate(queries):
        q_batch[i % procs].append(q)
    queries = None
    logging.info(f'  No. of batches: {len(q_batch))}')
    logging.info(f'  Queries per batch: {len(q_batch[0])}')
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

def write_table(idx: dict, outdir: str, qtax: str) -> None:
    """
    Write tab-delim table of taxonomy mappings to STDOUT.
    Params:
      idx: dict of {query_id: [taxname, count]}
      outdir: output directory
      qtax: either 'ncbi_taxonomy' or 'gtdb_taxonomy'
    Return:
      None
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, 'taxonomy_map_summary.tsv')
    ttax = 'ncbi_taxonomy' if qtax == 'gtdb_taxonomy' else 'gtdb_taxonomy'
    with open(outfile, 'w') as outF:
        outF.write('\t'.join([qtax, ttax, 'lca_frac',
                              'target_tax_level', 'lineage']) + '\n')
        for x in idx:
            for k,v in x.items():
                outF.write('\t'.join([k] + v) + '\n')
    logging.info(f'File written: {outfile}')

def rename(tax_idx: list, tax_queries: str, outdir: str, 
           column: int=1, header: bool=False) -> None:
    """
    Rename queries with new taxonomic classifications.
    All entries that cannot be re-named will be excluded in the output.
    Params:
      tax_idx: list of dicts of {query_id: [taxname, count]}
      tax_queries: text file of taxonomic queries
      outdir: output directory
      column: the column in the tax_queries file that contains the taxonomic classifications
      header: does the tax_queries file contain a header?
    Return:
      None
    """
    # Convert tax_idx to a simple index
    idx = {}  # {old_tax : new_tax}
    for x in tax_idx:
        for k,v in x.items():
            idx[k] = v[0]
    # Rename queries
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, 'queries_renamed.tsv')
    status = {'renamed' : 0, 'excluded' : 0}
    with gtdb2td.Utils.Open(tax_queries) as inF, open(outfile, 'w') as outF:
        for i,line in enumerate(inF):
            try:
                line = gtdb2td.Utils.Decode(line)
            except AttributeError:
                pass
            line = line.rstrip().split('\t')            
            if header is True and i == 0:
                pass
            else:
                try:
                    line[column-1] = idx[line[column-1]]
                    status['renamed'] += 1
                except KeyError:
                    status['excluded'] += 1
                    continue
                if line[column-1].lower() == 'unclassified':
                    status['renamed'] -= 1
                    status['excluded'] += 1
                    continue
            outF.write('\t'.join(line) + '\n')
    # status
    logging.info(f'File written: {outfile}')
    logging.info(f'  No. of queries renamed: {status["renamed"]}')
    logging.info(f'  No. of queries excluded: {status["excluded"]}')
            
def main(args: dict) -> None:
    """
    Main interface
    """
    # Load ncbi dmp files, if provided
    if args.names_dmp is not None and args.nodes_dmp is not None:
        ncbi_tax = gtdb2td.Dmp.load_dmp(args.names_dmp, args.nodes_dmp)
    else:
        ncbi_tax = None    
    # Load the metadata as graphs
    G = {'ncbi_taxonomy' : DiGraph_w_root(),
         'gtdb_taxonomy' : DiGraph_w_root()}
    for F in args.gtdb_metadata:
       load_gtdb_metadata(F, G, args.completeness, args.contamination,
                          no_prefix=args.no_prefix)
    # Query taxonomy
    idx = query_tax(args.tax_queries, G,
                    tax=args.query_taxonomy,
                    lca_frac = args.fraction,
                    max_tips = args.max_tips,
                    column = args.column,
                    header = args.header,
                    tax_graph = ncbi_tax,
                    prefix = args.prefix,
                    procs = args.procs, 
                    verbose = args.verbose)
    # Rename taxa
    if args.rename:
        rename(idx, args.tax_queries, args.outdir,
               column = args.column, header = args.header)
    # Write results
    write_table(idx, args.outdir, qtax=args.query_taxonomy)
             
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

