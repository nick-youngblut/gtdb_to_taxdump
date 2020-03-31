[![Travis-CI Build Status](https://travis-ci.org/nick-youngblut/gtdb_to_taxdump.svg?branch=master)](https://travis-ci.org/nick-youngblut/gtdb_to_taxdump)

gtdb_to_taxdump
===============

Convert GTDB taxonomy to NCBI taxdump format.

* Version: 0.1.2
* Authors:
  * Nick Youngblut <nyoungb2@gmail.com>
* Maintainers:
  * Nick Youngblut <nyoungb2@gmail.com>


# Summary

Convert GTDB taxonomy to NCBI taxdump format in order to
use the GTDB taxonomy with software that requires a
taxonomy in the taxdump format (eg., kraken2 or TaxonKit).

Note that the taxIDs are arbitrarily assigned and don't
match anything in the NCBI! Running `gtdb_to_taxdump` on
a different list of taxonomies (e.g., a different GTDB release)
will create different taxIDs.

# Citation 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3696964.svg)](https://doi.org/10.5281/zenodo.3696964)

# Install

No dependencies besides python >= 3.6

# Usage

See `gtdb_to_taxdump.py -h`

Example (GTDB release89):

```
./gtdb_to_taxdump.py \
  https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/ar122_taxonomy_r89.tsv \
  https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_taxonomy_r89.tsv \
  > taxID_info.tsv
```

You can add the taxIDs to a GTDB metadata table via the `--table` param. For example:

```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/ar122_metadata_r89.tsv
./gtdb_to_taxdump.py \
  --table ar122_metadata_r89.tsv \
  https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/ar122_taxonomy_r89.tsv \
  https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_taxonomy_r89.tsv \
  > taxID_info.tsv
```

# Extras

*  `gtdb_to_diamond.py`
  * Use to create a diamond database (with taxonomy) of genes from GTDB genomes
  * This can be used to taxonomically classify reads and amino acid sequences with diamond & GTDB

# GTDB website

https://data.ace.uq.edu.au/public/gtdb/data/releases/
