[![Travis-CI Build Status](https://travis-ci.org/nyoungblut/gtdb_to_taxdump.svg?branch=master)](https://travis-ci.org/nyoungblut/gtdb_to_taxdump)

gtdb_to_taxdump
===============

Convert GTDB taxonomy to NCBI taxdump format

* Version: 0.1.0
* Authors:
  * Nick Youngblut <nyoungb2@gmail.com>
* Maintainers:
  * Nick Youngblut <nyoungb2@gmail.com>


# Summary

Convert GTDB taxonomy to NCBI taxdump format in order to
use the GTDB taxonomy with software that requires a
taxonomy in the taxdump format (eg., kraken2 or TaxonKit).

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


# GTDB website

https://data.ace.uq.edu.au/public/gtdb/data/releases/
