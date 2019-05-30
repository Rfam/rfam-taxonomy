# Rfam Taxonomy

:warning: Work in progress, do not use in production :warning:

This repository contains the code and data for analysing the taxonomic distribution
of [Rfam](http://rfam.org) families.

The code uses the Rfam [public MySQL database](https://rfam.readthedocs.io/en/latest/database.html)
to compare the taxonomic domains of sequences from [seed alignments](https://rfam.readthedocs.io/en/latest/glossary.html)
and [full region](https://rfam.readthedocs.io/en/latest/glossary.html) hits.

:open_file_folder: The results are organised in several files in the [domains folder](./domains).
Each file contains seven columns:

1. `Family` = Rfam accession (e.g. RF00001)
2. `Domain` = Taxonomic domain where the family is found
3. `Seed domains` = All taxonomic domains from the seed alignment
4. `Full region domains` = All taxonomic domains from all the matches
5. `Rfam ID` = Rfam identifier (e.g. 5S_rRNA)
6. `Description` = Family description
7. `RNA type` = One of Rfam [RNA types](https://rfam.readthedocs.io/en/latest/searching-rfam.html#search-by-entry-type).

`Domain` can be:
- single domain (for example, Bacteria or Eukaryota) if the majority of hits (>=90%) are from the same domain both in seed and full region hits;
- `<seed domain>/<full region domain>` - if seed and full region domains are not the same, then both are listed. For example, `Viruses/Eukaryota` when seed alignment contains mostly Viruses and the full region hits contain mostly Eukaryotes);
- `Mixed` - if there is no single domain where the family occurs. For example, 5S rRNA RF00001 is expected to be found in Bacteria, Archaea, and Eukaryota.

:white_check_mark: View [summary](./domains/Readme.md) with the number of families observed in each domain.

The analysis is based on **Rfam 14.1**.

## Retrieving the data

The files can be retrieved directly from GitHub:

- [all-domains.csv](https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/all-domains.csv)
- [bacteria.csv](https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.csv)
- [archaea.csv](https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/archaea.csv)

 It is also possible to regenerate the files locally (see the **Installation** section below).

## Example use cases

If you are interested only in Rfam families that match Bacteria, you can use the [bacteria.csv](./domains/bacteria.csv) file. If you would like to do further processing, you can filter the [all-domains.csv](./domains/all-domains.csv) where the second column is `Bacteria`. For example, to eliminate any families that find hits outside Bacteria, make sure that the third and the fourth columns contain `Bacteria (100.0%)`.

--------------------------------------------------------------------------------

## Installation

### Requirements

- [virtualenv](https://virtualenv.pypa.io/en/latest/)
- [pip](https://pypi.org/project/pip/)

Clone or download this repository and run the following commands:

```
virtualenv ENV
source ENV/bin/activate
pip install -r requirements.txt
```

## Usage

```
# when running for the first time:
python rfam-taxonomy.py --precompute

# after precompute is done, run:
python rfam-taxonomy.py

# to see additional options:
python rfam-taxonomy.py --help
```
