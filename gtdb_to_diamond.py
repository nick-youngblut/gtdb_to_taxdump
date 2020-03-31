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
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdb_r89_rep_genomes.faa.tar.gz

Example "diamond makedb" run with output files:
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
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def copy_nodes(infile, outdir):
    """Simple copy of nodes.dmp file into the output directory
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
    """ Reading names.dmp file
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
    """ Getting .faa.gz files from the tarball
    """
    for tarinfo in members:
        x = os.path.splitext(tarinfo.name)
        if x[1] == '.gz':
            if os.path.splitext(x[0])[1] == '.faa':
                yield tarinfo

def faa_gz_index(directory='.', extension='.faa.gz'):
    """ Creating {accession:faa_file} index
    """
    regex = re.compile(r'_protein.faa.gz$')
    regexp = re.compile(r'^[^_]+_|_')    
    found = {}
    for dirpath, dirnames, files in os.walk(directory):
        for name in files:
            if name.lower().endswith(extension):
                accession = regexp.sub('', regex.sub('', name))
                found[accession] = os.path.join(dirpath, name)
    return found
                
def uncomp_tarball(tarball_file, tmp_dir):
    """ Extracting info from the tarball
    """
    # tmp dir
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    # extracting tarball
    logging.info('Extracting tarball: {}'.format(tarball_file))
    tar = tarfile.open(tarball_file)
    tar.extractall(path=tmp_dir, members=faa_gz_files(tar))
    tar.close()
    # listing files
    regex = re.compile(r'_protein.faa.gz$')
    faa_files = faa_gz_index(tmp_dir, '.faa.gz')
    msg = '  No. of .faa.gz files: {}'
    logging.info(msg.format(len(faa_files.keys())))
    return faa_files

def accession2taxid(names_dmp, faa_files, outdir):
    """ Creating accession2taxid table
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

def faa_merge(faa_files, outdir, gzip_out=False):
    """ Reading in, formatting, and merging all faa files
    """
    outfile = os.path.join(outdir, 'gtdb_all.faa')
    if gzip_out:
        _open = lambda x: gzip.open(x, 'wb')
        outfile += '.gz'
    else:
        _open = lambda x: open(x, 'w')
    logging.info('Formating & merging faa files...')
    with _open(outfile) as outF:
        for acc,faa_file in faa_files.items():            
            with gzip.open(faa_file) as inF:
                for line in inF:
                    line = line.decode('utf8')
                    if line.startswith('>'):
                        line = '>' + acc + ' ' + line.lstrip('>')
                    if gzip_out:
                        line = line.encode('utf-8')
                    outF.write(line)
    logging.info('  File written: {}'.format(outfile))
    
def main(args):
    """ Main interface
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
    if os.path.isdir(args.tmpdir):
        shutil.rmtree(args.tmpdir)
        logging.info('Temp-dir removed: {}'.format(args.tmpdir))
                   
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

