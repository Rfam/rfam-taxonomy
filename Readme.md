# Rfam Taxonomy

This repository contains the code and data for analysing the taxonomic distribution
of the [Rfam](http://rfam.org) families. The goal is to identify domain-specific
subsets of Rfam covariance models for annotating bacterial, eukaryotic,
and other genomes with the [Infernal](http://eddylab.org/infernal) software.

The code uses the Rfam [public MySQL database](https://rfam.readthedocs.io/en/latest/database.html)
to compare the taxonomic domains of sequences from the manually curated
[seed alignments](https://rfam.readthedocs.io/en/latest/glossary.html)
and the automatically identified [full region](https://rfam.readthedocs.io/en/latest/glossary.html) hits.

:open_file_folder: The results are organised in several files in the [domains folder](./domains).
Each file contains seven columns:

1. `Family` = Rfam accession (e.g. RF00001)
2. `Domain` = Taxonomic domain where the family is found
3. `Seed domains` = All taxonomic domains from the seed alignment
4. `Full region domains` = All taxonomic domains from full region hits
5. `Rfam ID` = Rfam identifier (e.g. 5S_rRNA)
6. `Description` = Family description
7. `RNA type` = One of Rfam [RNA types](https://rfam.readthedocs.io/en/latest/searching-rfam.html#search-by-entry-type).

`Domain` can be:
- a single domain (for example, _Bacteria_ or _Eukaryota_) if the majority of hits (>=90%) are from the same domain both in seed and full region hits;
- `<seed domain>/<full region domain>` - if seed and full region domains are not the same, then both are listed. For example, _Viruses/Eukaryota_ means that the seed alignment contains mostly Viruses and the full region hits contain mostly Eukaryotes);
- `Mixed` - if there is no single domain where the family occurs. For example, 5S rRNA [RF00001](http://rfam.org/family/RF00001) is expected to be found in _Bacteria_, _Archaea_, and _Eukaryota_.

:white_check_mark: View [summary](./domains/Readme.md) with the number of families observed in each domain.

The analysis is based on **Rfam 14.1**.

## Retrieving the data

The latest version of the files can be retrieved directly from GitHub using the following URL format:

- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/all-domains.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/archaea.csv

 It is also possible to download the data and use it locally or regenerate the files (see the **Installation** section below).

## Example use cases

- If you are interested in a subset of Rfam families that match Bacteria, you can use the [bacteria.csv](./domains/bacteria.csv) file. For example, the following command generates a `bacteria.cm` file with a subset of Rfam covariance models that can be used with the Infernal _cmscan_ program:

    ```
    curl https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.csv | \
    cut -f 1,1 -d ',' | \
    tail -n +2 | \
    cmfetch -o bacteria.cm -f Rfam.cm.gz -
    ```

    where _cmfetch_ is part of the Infernal suite and `Rfam.cm.gz` can be downloaded from `ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz`.

- You can also further process the [all-domains.csv](./domains/all-domains.csv) file. For example, to eliminate any families that find hits outside Bacteria, you can focus on rows where the second column is `Bacteria` and the third and the fourth columns contain `Bacteria (100.0%)`. Note that such a subset would ignore many important RNA families that detect some contamination in eukaryotic sequences.

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
python rfam-taxonomy.py --precompute-seed
python rfam-taxonomy.py --precompute-full

# after precompute is done, run:
python rfam-taxonomy.py

# to see additional options:
python rfam-taxonomy.py --help
```

## Feedback

Feel free to create GitHub issues to ask questions or provide feedback.
Pull requests are also welcome.
