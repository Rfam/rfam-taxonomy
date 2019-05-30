# Rfam Taxonomy

:warning: Work in progress, do not use in production :warning:

The file [domains.csv](./domains.csv) contains seven columns:

- `Family` = Rfam accession (e.g. RF00001)
- `Domain` = Taxonomic domain where the family is found
- `Seed domains` = All taxonomic domains from the seed alignment
- `Full region domains` = All taxonomic domains from all the matches
- `Rfam ID` = Rfam identifier (e.g. 5S_rRNA)
- `Description` = Family description
- `RNA type` = One of Rfam [RNA types](https://rfam.readthedocs.io/en/latest/searching-rfam.html#search-by-entry-type)

Example:

```
Family,Domain,Seed domains,Full region domains
RF00001,Mixed,"Bacteria (48.6%), Eukaryota (45.51%), Archaea (5.9%)","Eukaryota (87.59%), Bacteria (12.0%), Archaea (0.4%)"
RF00002,Eukaryota,Eukaryota (100.0%),"Eukaryota (99.54%), Bacteria (0.46%)"
RF00003,Eukaryota,Eukaryota (100.0%),Eukaryota (100.0%)
RF00004,Eukaryota,"Eukaryota (99.04%), unclassified sequences (0.96%)","Eukaryota (99.98%), Viruses (0.01%), unclassified sequences (0.01%)"
RF00005,Mixed,"Eukaryota (78.2%), Bacteria (17.3%), Archaea (3.46%), Viruses (1.05%)","Eukaryota (73.83%), Bacteria (25.02%), Archaea (0.95%), Viruses (0.2%)"
```

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
