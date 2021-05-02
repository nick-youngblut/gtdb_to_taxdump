![gtdb_to_taxdump](https://github.com/nick-youngblut/gtdb_to_taxdump/workflows/gtdb_to_taxdump/badge.svg)

gtdb_to_taxdump
===============

Convert GTDB taxonomy to NCBI taxdump format.

* Version: 0.1.6
* Authors:
  * Nick Youngblut <nyoungb2@gmail.com>
* Maintainers:
  * Nick Youngblut <nyoungb2@gmail.com>

# WARNING

> There was a serious bug with `ncbi-gtdb_map.py` prior to version 0.1.5.
  Many of the taxonomic classifications are likely incorrect.
  Please re-run the analysis. I'm sorry for any inconvenience.

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

* `python >= 3.6`
* Only for `gtdb_to_diamond.py`
  * `diamond`
* For `ncbi-gtdb_map.py` & `lineage2taxid.py`
  * `networkx >= 2.4`

# Usage

See `gtdb_to_taxdump.py -h`

Example (GTDB release95):

```
./gtdb_to_taxdump.py \
  https://data.gtdb.ecogenomic.org/releases/release202/202.0/ar122_taxonomy_r202.tsv.gz \
  https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_taxonomy_r202.tsv.gz \
  > taxID_info.tsv
```

Example (GTDB release95):

```
./gtdb_to_taxdump.py \
  https://data.gtdb.ecogenomic.org/releases/release95/95.0/ar122_taxonomy_r95.tsv.gz \
  https://data.gtdb.ecogenomic.org/releases/release95/95.0/bac120_taxonomy_r95.tsv.gz \
  > taxID_info.tsv
```

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

* `ncbi-gtdb_map.py`
  * Map taxonomic classifications between the NCBI and GTDB taxonomies
  * Taxonomy mapping is done via the GTDB metadata table (includes both NCBI & GTDB taxonomies)
  * Mapping can go either way: `NCBI => GTDB` or `GTDB => NCBI` (see `--query-taxonomy`)
  * The user can select "fuzzy" mappings, in which the mapping isn't a perfect 1-to-1
  * See the script help docs for examples on usage
* `gtdb_to_diamond.py`
  * Use to create a diamond database (with taxonomy) of genes from GTDB genomes
  * This can be used to taxonomically classify reads and amino acid sequences with diamond & GTDB
  * See the script help docs for examples on usage
* `lineage2taxid.py`
  * Use to get the taxids for a set of lineages (eg., from GTDB-Tk) 
  * This is somewhat of a reverse of `gtdb_to_taxdump.py`
  * The script can work with GTDB or NCBI tax dumps
* `./uniref_utils/`
  * `unirefxml2clust50-90idx.py`
    * Mapping UniRef90 cluster IDs to UniRef50 cluster IDs
  * `unirefxml2fasta.py`
    * Creating gene sequence fasta from xml
  * `unirefxml2tax.py`
    * Converting the xml to a table of gene taxonomies

# GTDB website

https://data.ace.uq.edu.au/public/gtdb/data/releases/

