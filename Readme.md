# Rfam Taxonomy

:warning: Work in progress, do not use in production :warning:

The file [domains.csv](./domains.csv) contains seven columns:

1. `Family` = Rfam accession (e.g. RF00001)
2. `Domain` = Taxonomic domain where the family is found
3. `Seed domains` = All taxonomic domains from the seed alignment
4. `Full region domains` = All taxonomic domains from all the matches
5. `Rfam ID` = Rfam identifier (e.g. 5S_rRNA)
6. `Description` = Family description
7. `RNA type` = One of Rfam [RNA types](https://rfam.readthedocs.io/en/latest/searching-rfam.html#search-by-entry-type)

`Domain` can be:
- single domain (like Bacteria or Eukaryota) if the majority of hits (>=90%) are from the same domain both in seed and full region hits;
- `<domain>/<domain>` - if seed and full region domains are not the same, both are listed (`Viruses/Eukaryota` when seed alignment contains mostly Viruses and the full region hits contain mostly Eukaryotes);
- `Mixed` if there is no dominant domain

:white_check_mark: View [summary](./domains/summary.md) with the number of families observed in each domain.

## Example use cases

If you are interested only in families that match Bacteria, you can filter the file `domains.csv` where the second column is `Bacteria`. To further eliminate any families that find hits outside Bacteria, make sure the fourth column is `Bacteria (100.0%)`.


Please note that some families are expected to be found in multiple domains (for example RF00001).

For example, [RF01390](http://rfam.org/family/RF01390) has only Bacteria in the seed but 91% of full region hits are Eukaryotes. In this example this could be due to contamination as this RNA can be expressed [in pathogenic bacteria located inside eukaryotic cells](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4228915/).

The analysis is based on **Rfam 14.1**.

--------------------------------------------------------------------------------

## Installation

```
virtualenv ENV
source ENV/bin/activate
pip install -r requirements.txt
```

## Usage

```
python rfam-taxonomy.py

# to see additional options:
python rfam-taxonomy.py --help
```
