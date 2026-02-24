# Rfam Taxonomy

Based on **Rfam 15.1** (December 2024). See [releases](https://github.com/Rfam/rfam-taxonomy/releases) for previous versions.

This repository contains the code and data for analysing the taxonomic distribution
of the [Rfam](https://rfam.org) families. The goal is to identify taxonomic group-specific
subsets of Rfam covariance models for annotating bacterial, eukaryotic,
and other genomes or sequence datasets with the [Infernal](http://eddylab.org/infernal) software.

The code uses the Rfam [public MySQL database](https://rfam.readthedocs.io/en/latest/database.html)
to compare the taxonomic groups (domains and subgroups) of sequences from the manually curated
[seed alignments](https://rfam.readthedocs.io/en/latest/glossary.html)
and the automatically identified [full region](https://rfam.readthedocs.io/en/latest/glossary.html) hits.

:open_file_folder: The results are organised in several files in the [domains folder](./domains).
Each group (domain or subgroup) has both a `.csv` file and a `.clanin` file:

- **CSV files** contain seven columns with family and taxonomic information
- **Clanin files** list Rfam clan membership for families in that group, useful for `cmscan --clanin` when using modified CM files

CSV files contain seven columns:

1. `Family` = Rfam accession (e.g. RF00001)
2. `Domain` = Taxonomic classification where the family is found (:grey_exclamation: this is the most important column)
3. `Seed groups` = All taxonomic groups (domains and subgroups) from the seed alignment
4. `Full region groups` = All taxonomic groups (domains and subgroups) from full region hits
5. `Rfam ID` = Rfam identifier (e.g. 5S_rRNA)
6. `Description` = Family description
7. `RNA type` = One of Rfam [RNA types](https://rfam.readthedocs.io/en/latest/searching-rfam.html#search-by-entry-type).

## Rules for the `Domain` field

The `Domain` field (column 2) is determined by the following rules:

- **Single group** (e.g., _Bacteria_, _Eukaryota_): Used if the majority of hits (≥90%) are from the same group in both seed and full region hits.
- **A+B**: If multiple groups are each ≥90% within a single dataset (seed or full), they are joined with `+` (e.g., _Eukaryota+Fungi_). This occurs with parent/subgroup relationships where the subgroup's sequences contribute to both the parent and subgroup counts, allowing both to exceed 90%.
- **A/B**: If seed and full region groups are different, the field is `<seed group>/<full region group>`. For example, _Viruses/Eukaryota_ means the seed alignment is mostly Viruses and the full region hits are mostly Eukaryota. This also applies to parent/subgroup pairs when they differ between seed and full.
- **Mixed**: Used if there is no single group where the family occurs (i.e., no group is above the cutoff and the distribution is ambiguous). For example, 5S rRNA [RF00001](http://rfam.org/family/RF00001) is found in _Bacteria_, _Archaea_, and _Eukaryota_.
- **A/Mixed** or **Mixed/B**: If one of the seed or full region groups is `Mixed`, the field is shown as `<seed group>/Mixed` or `Mixed/<full region group>`.
- **A/No Data** or **No Data/B**: If there is no data for a family in either the seed or full region, the field is shown as `<seed group>/No Data`, `No Data/<full region group>`, or `No Data` if both lack data.

### Parent-child domain relationships

Some taxonomic groups are considered subgroups of major domains. As of now, parent domains must be one of the major taxonomic domains: _Archaea_, _Bacteria_, _Eukaryota_, _Viruses_, _Viroids_, or _Other_.

**Biological domains (top-level):**
- _Archaea_
- _Bacteria_
- _Eukaryota_
- _Viruses_
- _Viroids_
- _unclassified sequences_
- _Other_

**Currently defined subgroups:**

- **Eukaryota subgroups:** _Fungi_, _Viridiplantae_, _Metazoa_, _Amoebozoa_, _Alveolata_, _Stramenopiles_
- **Bacteria subgroups:** _Proteobacteria_, _Firmicutes_, _Actinobacteria_, _Bacteroidetes_, _Cyanobacteria_
- **Archaea subgroups:** _Euryarchaeota_, _Crenarchaeota_, _Thaumarchaeota_
- **Virus subgroups:** _Flaviviridae_, _Coronaviridae_, _Retroviridae_, _Herpesviridae_, _Picornaviridae_, _Orthomyxoviridae_

When both a parent domain and its subgroup each account for ≥90% of the sequences in both the seed and full datasets (because subgroup sequences are counted toward both the subgroup and its parent), and the resulting major group labels for seed and full are the same and use this pattern, the `Domain` field is shown as `Parent+Subgroup` (e.g., _Eukaryota+Fungi_). This helps clarify cases where a family is strongly represented in both a major domain and a specific subgroup.

In the `Seed groups` and `Full region groups` columns (columns 3 and 4), subgroups are displayed with their parent prefix using a colon separator to indicate the hierarchical relationship (e.g., _Eukaryota:Fungi_). The percentages for subgroups are also included in their parent domain's percentage.

:white_check_mark: View [summary](./domains/Readme.md) with the number of families observed in each group.

## Retrieving the data

The latest version of the files can be retrieved directly from GitHub using the following URL format:

- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/all-domains.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/archaea.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/archaea.clanin
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.clanin
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/eukaryota.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/eukaryota.clanin
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/viruses.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/viruses.clanin
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/viroids.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/viroids.clanin
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/unclassified-sequences.csv
- https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/unclassified-sequences.clanin

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

    The same approach works for any group (e.g. archaea or fungi) by just replacing `bacteria` in the command above with the name of the group.

    **Note for subgroups:** Subgroup files like [fungi.csv](./domains/fungi.csv) contain two types of families: (1) families that are ≥90% from that subgroup (Domain field = `Fungi`), and (2) families that are primarily from the parent domain but have significant (>5%) representation from the subgroup (Domain field = `Eukaryota` or containing the parent name). To extract only families that are ≥90% from a subgroup, filter the Domain field (column 2):

    ```
    # Get only families that are ≥90% Fungi
    curl https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/fungi.csv | \
    awk -F',' '$2 ~ /(^|\+)Fungi(\+|$)/ {print $1}' | \
    tail -n +2 | \
    cmfetch -o fungi-90pct.cm -f Rfam.cm.gz -
    ```

- When using group-specific CM files (created as shown above), you should use the corresponding `.clanin` file with `cmscan --clanin` to ensure proper clan competition filtering. The `--clanin` option must be used with `--fmt 2` and `--tblout`. For example:

    ```
    cmscan --cpu 4 --fmt 2 --clanin bacteria.clanin --tblout results.tbl bacteria.cm sequences.fa
    ```

    The `.clanin` files contain clan membership information for families present in each group, ensuring that `cmscan` correctly handles competing models from the same clan when using modified CM subsets.

    **Additional clan-related options:**
    - `--oclan`: Only mark overlaps between models in the same clan (requires `--clanin`, `--fmt 2`, and `--tblout`)
    - `--oskip`: Omit lower-scoring overlapping hits from the tabular output file (requires `--fmt 2` and `--tblout`). When combined with `--oclan`, only omits overlaps within the same clan.

    Example using all clan options:
    ```
    cmscan --cpu 4 --fmt 2 --clanin bacteria.clanin --oclan --oskip --tblout results.tbl bacteria.cm sequences.fa
    ```

- You can also further process the [all-domains.csv](./domains/all-domains.csv) file. For example, to eliminate any families that find hits outside Bacteria, you can focus on rows where the second column is `Bacteria` and the third and the fourth columns (the group percentages) contain `Bacteria (100.0%)`. Note that such a subset would ignore many important RNA families that detect some contamination in eukaryotic sequences.

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

    # this will generate both .csv and .clanin files for each group
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
