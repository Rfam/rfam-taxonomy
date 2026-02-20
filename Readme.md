# Rfam Taxonomy

Based on **Rfam 15.1** (December 2024). See [releases](https://github.com/Rfam/rfam-taxonomy/releases) for previous versions.

This repository contains the code and data for analysing the taxonomic distribution
of the [Rfam](https://rfam.org) families. The goal is to identify domain-specific
subsets of Rfam covariance models for annotating bacterial, eukaryotic,
and other genomes with the [Infernal](http://eddylab.org/infernal) software.

The code uses the Rfam [public MySQL database](https://rfam.readthedocs.io/en/latest/database.html)
to compare the taxonomic domains of sequences from the manually curated
[seed alignments](https://rfam.readthedocs.io/en/latest/glossary.html)
and the automatically identified [full region](https://rfam.readthedocs.io/en/latest/glossary.html) hits.

:open_file_folder: The results are organised in several files in the [domains folder](./domains).
Each domain has both a `.csv` file and a `.clanin` file:

- **CSV files** contain seven columns with family and taxonomic information
- **Clanin files** list Rfam clan membership for families in that domain, useful for `cmscan --clanin` when using modified CM files

CSV files contain seven columns:

1. `Family` = Rfam accession (e.g. RF00001)
2. `Domain` = Taxonomic domain where the family is found (:grey_exclamation: this is the most important column)
3. `Seed domains` = All taxonomic domains from the seed alignment
4. `Full region domains` = All taxonomic domains from full region hits
5. `Rfam ID` = Rfam identifier (e.g. 5S_rRNA)
6. `Description` = Family description
7. `RNA type` = One of Rfam [RNA types](https://rfam.readthedocs.io/en/latest/searching-rfam.html#search-by-entry-type).

`Domain` can be:
- a single domain (for example, _Bacteria_ or _Eukaryota_) if the majority of hits (>=90%) are from the same domain both in seed and full region hits;
- `<seed domain>/<full region domain>` - if seed and full region domains are not the same, then both are listed. For example, _Viruses/Eukaryota_ means that the seed alignment contains mostly Viruses and the full region hits contain mostly Eukaryotes);
- `Mixed` - if there is no single domain where the family occurs. For example, 5S rRNA [RF00001](http://rfam.org/family/RF00001) is expected to be found in _Bacteria_, _Archaea_, and _Eukaryota_.
- `<seed region domain>/Mixed` or `Mixed/<full region domain>` - For example, Bacterial SSU [RF00177](http://rfam.org/family/RF00177) has only Bacteria in the seed alignment but the full region hits also contain Eukaryota because the mitochondrial and plastid SSU is similar to the bacterial SSU and is expected to match the bacterial model.

:white_check_mark: View [summary](./domains/Readme.md) with the number of families observed in each domain.

## Retrieving the data

The latest version of the files can be retrieved directly from GitHub using the following URL format:

- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/all-domains.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.clanin
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/archaea.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/archaea.clanin
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/viruses.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/viruses.clanin

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

- When using domain-specific CM files (created as shown above), you should use the corresponding `.clanin` file with `cmscan --clanin` to ensure proper clan competition filtering. The `--clanin` option must be used with `--fmt 2` and `--tblout`. For example:

    ```
    cmscan --cpu 4 --fmt 2 --clanin bacteria.clanin --tblout results.tbl bacteria.cm sequences.fa
    ```

    The `.clanin` files contain clan membership information for families present in each domain, ensuring that `cmscan` correctly handles competing models from the same clan when using modified CM subsets.

    **Additional clan-related options:**
    - `--oclan`: Only mark overlaps between models in the same clan (requires `--clanin`, `--fmt 2`, and `--tblout`)
    - `--oskip`: Omit lower-scoring overlapping hits from the tabular output file (requires `--fmt 2` and `--tblout`). When combined with `--oclan`, only omits overlaps within the same clan.

    Example using all clan options:
    ```
    cmscan --cpu 4 --fmt 2 --clanin bacteria.clanin --oclan --oskip --tblout results.tbl bacteria.cm sequences.fa
    ```

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

## Updating the data

After each Rfam release, the data in this repo need to be updated locally and
pushed to GitHub.

1. Generate new data

    ```
    # when running for the first time (needs to run in this order):
    python rfam-taxonomy.py --precompute-full
    python rfam-taxonomy.py --precompute-seed

    # after precompute is done, run:
    python rfam-taxonomy.py

    # this will generate both .csv and .clanin files for each domain
    # clan information is retrieved directly from the Rfam database

    # to see additional options:
    python rfam-taxonomy.py --help
    ```

1. Review the changes

    The results must be manually reviewed before committing the new files by checking the difference between the old and the new versions using git.

    It is normal for the values in the 3rd and 4th columns to change but `Domain`, the 2nd column, should stay stable unless the affected family has been significantly updated.

1. Update release info in Readme

1. Create new GitHub release

## Feedback

Feel free to create [GitHub issues](https://github.com/Rfam/rfam-taxonomy/issues) to ask questions or provide feedback.
Pull requests are also welcome.
