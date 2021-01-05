#!/usr/bin/env python
from __future__ import print_function
# batteries
import os
import re
import sys
import gzip
import glob
import shutil
import argparse
import logging
import urllib.request
import codecs
import tarfile
from collections import OrderedDict


desc = 'Converting GTDB taxonomy to input for "diamond makedb --taxonmap"'
epi = """DESCRIPTION:
Convert Genome Taxonomy Database (GTDB) representative genome
gene amino acid sequences to the input files required for 
"diamond makedb --taxonmap", with allows for taxonomic identification
with diamond (LCA-based method).

Example of getting a GTDB faa fasta tarball:
  
  wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/genomic_files_reps/gtdb_proteins_aa_reps_r95.tar.gz

Example extraction & formatting of faa files from tarball:

  gtdb_to_diamond.py gtdb_proteins_aa_reps_r95.tar.gz names.dmp nodes.dmp

Example "diamond makedb" run with gtdb_to_diamond.py output files:

  diamond makedb --in $OUTDIR/gtdb_all.faa --db gtdb.dmnd --taxonmap $OUTDIR/accession2taxid.tsv --taxonnodes $OUTDIR/nodes.dmp --taxonnames $OUTDIR/names.dmp
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('faa_tarball', metavar='faa_tarball', type=str,
                    help='tarball of GTDB ref genome gene animo acid data files')
parser.add_argument('names_dmp', metavar='names_dmp', type=str,
                    help='taxdump names.dmp file (eg., from gtdb_to_taxdump.py)')
parser.add_argument('nodes_dmp', metavar='nodes_dmp', type=str,
                    help='taxdump nodes.dmp file (eg., from gtdb_to_taxdump.py)')
parser.add_argument('-o', '--outdir', type=str, default='gtdb_to_diamond',
                    help='Output directory (Default: %(default)s)')
parser.add_argument('-t', '--tmpdir', type=str, default='gtdb_to_diamond_TMP',
                    help='Temporary directory (Default: %(default)s)')
parser.add_argument('-g', '--gzip', action='store_true', default=False,
                    help='gzip output fasta? (Default: %(default)s)')
parser.add_argument('-k', '--keep-temp', action='store_true', default=False,
                    help='Keep temporary output? (Default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def copy_nodes(infile, outdir):
    """
    Simple copy of nodes.dmp file into the output directory
    """
    logging.info('Read nodes.dmp file: {}'.format(infile))
    outfile = os.path.join(outdir, 'nodes.dmp')
    if infile == outfile:
        raise IOError('Input == Output: {} <=> {}'.format(infile, outfile))
    with open(infile) as inF, open(outfile, 'w') as outF:
        for line in inF:
            outF.write(line)
    logging.info('File written: {}'.format(outfile))

def read_names_dmp(infile, outdir):
    """ 
    Reading names.dmp file
    """
    outfile = os.path.join(outdir, 'names.dmp')
    regex = re.compile(r'\t\|\t')
    regexp = re.compile(r'^[^_]+_|_')
    names_dmp = {}
    logging.info('Reading dumpfile: {}'.format(infile))
    with open(infile) as inF, open(outfile, 'w') as outF:
        for line in inF:
            line = regex.split(line.rstrip())
            if len(line) >= 2:
                line[1] = regexp.sub('', line[1])  # accession
                names_dmp[line[1]] = line[0]
            outF.write('\t|\t'.join(line) + '\n')
            
    logging.info('  File written: {}'.format(outfile))
    msg = '  No. of accession<=>taxID pairs: {}'
    logging.info(msg.format(len(names_dmp.keys())))
    return names_dmp

def faa_gz_files(members):
    """ 
    Getting .faa.gz files from the tarball
    """
    for tarinfo in members:
        for ext in ('.faa.gz', '.faa'):
            if tarinfo.name.endswith(ext):
                yield tarinfo

def faa_gz_index(directory='.', extensions=['.faa', '.faa.gz']):
    """ 
    Creating {accession:faa_file} index from extracted tarball files
    """
    extensions = set(extensions)
    regex = re.compile(r'(_protein\.faa\.gz|_protein\.faa)$')
    regexp = re.compile(r'^[^_]+_|_')    
    found = {}
    for dirpath, dirnames, files in os.walk(directory):
        for name in files:
            for ext in extensions:
                if name.lower().endswith(ext):
                    accession = regexp.sub('', regex.sub('', name))
                    found[accession] = os.path.join(dirpath, name)
                    continue
    return found
                
def uncomp_tarball(tarball_file, tmp_dir):
    """ 
    Extracting info from the tarball
    """
    # tmp dir
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    # extracting tarball
    logging.info('Extracting tarball: {}'.format(tarball_file))
    logging.info('  Extracting to: {}'.format(tmp_dir))
    tar = tarfile.open(tarball_file)
    tar.extractall(path=tmp_dir, members=faa_gz_files(tar))
    tar.close()
    # listing files
    faa_files = faa_gz_index(tmp_dir, ['.faa', '.faa.gz'])
    n_files = len(faa_files.keys())
    msg = '  No. of .faa(.gz) files: {}'
    logging.info(msg.format(n_files))
    if n_files == 0:
        logging.warning('  No .faa(.gz) files found!')
    return faa_files

def accession2taxid(names_dmp, faa_files, outdir):
    """ 
    Creating accession2taxid table
    """
    outfile = os.path.join(outdir, 'accession2taxid.tsv')    
    logging.info('Creating accession2taxid table...')
    with open(outfile, 'w') as outF:
        header = ['accession', 'accession.version', 'taxid', 'gi']
        outF.write('\t'.join(header) + '\n')
        for accession,faa_file in faa_files.items():
            try:
                taxID = names_dmp[accession]
            except KeyError:
                msg = 'Cannot find {} accession in names.dmp'
                raise KeyError(msg.format(accession))
            acc_base = os.path.splitext(accession)[0]
            line = [acc_base, accession, str(taxID), '']
            outF.write('\t'.join(line) + '\n')
    logging.info('  File written: {}'.format(outfile))

def _open(filename, io='r'):
    """
    Read/Write (gzip'ed) files
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, io + 'b')
    else:
        return open(filename, io)
    
def faa_merge(faa_files, outdir, gzip_out=False):
    """ 
    Reading in, formatting, and merging all faa files
    """
    outfile = os.path.join(outdir, 'gtdb_all.faa')
    if gzip_out:
        outfile += '.gz'
    
    logging.info('Formating & merging faa files...')
    seq_cnt = 0
    with _open(outfile, 'w') as outF:
        for acc,faa_file in faa_files.items():            
            with _open(faa_file, 'r') as inF:
                for line in inF:
                    try:
                        line = line.decode('utf8')
                    except AttributeError:
                        pass
                    if line.startswith('>'):
                        line = '>' + acc + ' ' + line.lstrip('>')
                        seq_cnt += 1
                    if gzip_out:
                        line = line.encode('utf-8')
                    outF.write(line)
    logging.info('  File written: {}'.format(outfile))
    logging.info('  No. of seqs. written: {}'.format(seq_cnt))
    
def main(args):
    """ 
    Main interface
    """
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    # copying nodes
    copy_nodes(args.nodes_dmp, args.outdir)
    # reading in names.dmp
    names_dmp = read_names_dmp(args.names_dmp, args.outdir)
    # uncompressing tarball of faa fasta files
    faa_files = uncomp_tarball(args.faa_tarball, args.tmpdir)
    # create accession2taxid
    accession2taxid(names_dmp, faa_files, args.outdir)
    # creating combined faa fasta
    faa_merge(faa_files, args.outdir, gzip_out=args.gzip)
    # clean up
    if not args.keep_temp and os.path.isdir(args.tmpdir):
       shutil.rmtree(args.tmpdir)
       logging.info('Temp-dir removed: {}'.format(args.tmpdir))
                   
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

