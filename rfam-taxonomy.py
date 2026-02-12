"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""


import os
import csv
import re

import click

from scripts.rfam_db import get_rfam_families, get_taxonomy_info


DATA_SEED_PATH = 'data-seed'
DATA_FULL_REGION_PATH = 'data-full-region'

DOMAINS = sorted([
    'Archaea',
    'Bacteria',
    'Eukaryota',
    'Fungi',
    'unclassified sequences',
    'Viruses',
    'Viroids',
    'Other',
])
DOMAIN_CUTOFF = 90 # at least 90% of sequences must be from this domain
FUNGI_THRESHOLD = 5  # include families with at least 5% Fungi in full regions

# Mapping of subgroups to their parent domains
# To add a new subgroup, add it here and to the DOMAINS list above
SUBGROUP_PARENT = {
    'Fungi': 'Eukaryota',
    # Add more subgroup-parent pairs here, e.g.:
    # 'Plants': 'Eukaryota',
    # 'Alphaproteobacteria': 'Bacteria',
}

WHITELIST = [
    'RF00001', # 5S
    'RF00005', # tRNA
    'RF01852', # tRNA-Sec
]


def precompute_taxonomic_information(analysis_type):
    """
    Store csv files for each family (lineage, count, NCBI taxid).

    Example:
    Bacteria; Acidobacteria; Acidobacteriales; Acidobacteriaceae; Acidobacterium.,1,240015
    Bacteria; Actinobacteria; Acidimicrobidae; Acidimicrobiales; Acidimicrobineae; Acidimicrobiaceae; Acidimicrobium.,1,525909
    Bacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Catenulisporineae; Catenulisporaceae; Catenulispora.,1,479433
    """
    print('Retrieving data from the public Rfam database')

    if analysis_type == 'seed':
        DATA_PATH = DATA_SEED_PATH
    elif analysis_type == 'full':
        DATA_PATH = DATA_FULL_REGION_PATH
    os.system('mkdir -p {}'.format(DATA_PATH))

    for family in get_rfam_families():
        print(family['rfam_acc'])
        with open(os.path.join(DATA_PATH, '{}.csv'.format(family['rfam_acc'])), 'w') as csvfile:
            for row in get_taxonomy_info(family['rfam_acc'], analysis_type):
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(row)
    print('Done')


def get_taxonomic_distribution(rfam_acc, DATA_PATH):
    """
    Calculate the percentage of hits from each domain for a family.

    Example:
    {'Eukaryota': 45.51, 'Bacteria': 48.6, 'Other': 0.0, 'Viruses': 0.0, 'unclassified sequences': 0.0, 'Viroids': 0.0, 'Archaea': 5.9}
    """
    data = {}
    for domain in DOMAINS:
        data[domain] = 0
    total = 0
    num_rows = 0
    with open(os.path.join(DATA_PATH, '{}.csv'.format(rfam_acc)), 'r') as f:
        csvreader = csv.reader(f)
        for row in csvreader:
            num_rows += 1
            tax_string, count, _ = row
            count = int(count)
            taxon = tax_string.split(';')[0]

            if '.' in taxon:
                taxon = taxon.replace('.', '') # example: Bacteria.
            if taxon == 'Unclassified':
                taxon = 'unclassified sequences'

            # Check if Fungi appears in the taxonomy lineage (sub-kingdom of Eukaryota)
            if 'Fungi' in tax_string:
                taxon = 'Fungi'

            if taxon in DOMAINS:
                data[taxon] += count
            else:
                data['Other'] += count
            total += count
    if num_rows == 0:
        return 'NO_DATA'
    if total != 0:
        for domain in DOMAINS:
            data[domain] = round(data[domain]*100.0/total, 2)
    return data


def get_major_domain(data, cutoff):
    """
    Find the prevalent domain (for example, Eukaryota):

    {'Eukaryota': 100.0, 'Other': 0.0, 'Viruses': 0.0, 'unclassified sequences': 0.0, 'Viroids': 0.0, 'Archaea': 0.0, 'Bacteria': 0.0}
    """
    if data == 'NO_DATA':
        return 'No Data'
    # Find all domains above cutoff
    major_domains = [domain for domain, value in data.items() if value >= cutoff]
    if len(major_domains) == 1:
        return major_domains[0]
    elif len(major_domains) > 1:
        return '+'.join(sorted(major_domains))
    else:
        found_domains = [domain for domain, value in data.items() if value > 0]
        # Guard against empty found_domains
        if not found_domains:
            return 'Mixed'
        # Generalized special case: only a parent and its subgroup present
        for subgroup, parent in SUBGROUP_PARENT.items():
            if set(found_domains) <= set([parent, subgroup]):
                # Assign to the one with the higher count
                if data[parent] > data[subgroup]:
                    return parent
                elif data[subgroup] > data[parent]:
                    return subgroup
                else:
                    # If equal, return both joined
                    return '+'.join(sorted([parent, subgroup]))
        # if only two domains and one of them is `unclassified`, consider the other one major domain
        if len(found_domains) == 2 and 'unclassified sequences' in found_domains:
            found_domains.remove('unclassified sequences')
            return found_domains.pop()
        else:
            return 'Mixed'


def get_domains(data):
    """
    List all domains in which a family has been observed.
    """
    if data == 'NO_DATA':
        return 'No Data'

    # Make a copy so we don't mutate the original
    data = data.copy()

    # For each subgroup, add its percentage to its parent
    for subgroup, parent in SUBGROUP_PARENT.items():
        if parent in data and subgroup in data:
            data[parent] += data[subgroup]
            # Cap at 100%
            if data[parent] > 100.0:
                data[parent] = 100.0

    output = []
    for domain, proportion in sorted(data.items(), key=lambda x: x[1], reverse=True):
        if proportion > 0:
            # Show subgroups with their parent prefix
            display_name = domain
            if domain in SUBGROUP_PARENT:
                display_name = '{}/{}'.format(SUBGROUP_PARENT[domain], domain)
            output.append('{} ({:.2f}%)'.format(display_name, proportion))
    return ', '.join(output)


def get_fungi_percentage(full_domains_str):
    """
    Extract the Fungi percentage from a full_domains string.

    Example input: 'Bacteria (0.25%), Eukaryota/Fungi (9.66%)'
    Returns: 9.66

    If Fungi is not present, returns 0.0
    """
    # Match "Fungi (x%)" possibly preceded by a parent prefix like "Eukaryota/Fungi"
    # within a single comma-separated domain entry.
    match = re.search(r'[^,]*\\bFungi \\(([0-9.]+)%\\)', full_domains_str)
    if match:
        return float(match.group(1))
    return 0.0


def analyse_seed_full_taxonomic_distribution(family, cutoff):
    """
    Compare domains observed in seed alignments and full region hits.
    """
    seed = get_taxonomic_distribution(family['rfam_acc'], DATA_SEED_PATH)
    full = get_taxonomic_distribution(family['rfam_acc'], DATA_FULL_REGION_PATH)

    major_domain_seed = get_major_domain(seed, cutoff)
    seed_domains = get_domains(seed)

    major_domain_full = get_major_domain(full, cutoff)
    full_domains = get_domains(full)

    # Rule 1: If both are the same and not Mixed/No Data, use that value
    if major_domain_seed and major_domain_seed == major_domain_full and major_domain_seed not in ('Mixed', 'No Data'):
        domain_field = major_domain_seed
    # Rule 2: If either is Mixed or No Data, use the appropriate combination
    elif major_domain_seed in ('Mixed', 'No Data') or major_domain_full in ('Mixed', 'No Data'):
        # Handle all combinations of Mixed and No Data
        if major_domain_seed == 'No Data' and major_domain_full == 'No Data':
            domain_field = 'No Data'
        elif major_domain_seed == 'Mixed' and major_domain_full == 'Mixed':
            domain_field = 'Mixed'
        elif major_domain_seed == 'No Data' and major_domain_full == 'Mixed':
            domain_field = 'Mixed'
        elif major_domain_seed == 'Mixed' and major_domain_full == 'No Data':
            domain_field = 'Mixed'
        elif major_domain_seed == 'No Data':
            domain_field = 'No Data/{}'.format(major_domain_full)
        elif major_domain_full == 'No Data':
            domain_field = '{}/No Data'.format(major_domain_seed)
        elif major_domain_seed == 'Mixed':
            domain_field = 'Mixed/{}'.format(major_domain_full)
        else:  # major_domain_full == 'Mixed'
            domain_field = '{}/Mixed'.format(major_domain_seed)
    else:
        # Rule 3: Handle A+B patterns and standard A/B format
        # If one is A+B and the other is A (or any component of A+B), use A+B
        if (isinstance(major_domain_seed, str) and '+' in major_domain_seed and major_domain_full in major_domain_seed.split('+')):
            domain_field = major_domain_seed
        elif (isinstance(major_domain_full, str) and '+' in major_domain_full and major_domain_seed in major_domain_full.split('+')):
            domain_field = major_domain_full
        # Otherwise use standard seed/full format
        else:
            domain_field = '{}/{}'.format(major_domain_seed, major_domain_full)
    return [
        family['rfam_acc'],
        domain_field,
        seed_domains,
        full_domains,
        family['rfam_id'],
        family['description'],
        family['type'],
    ]


def write_output_files(data):
    """
    Generate output files.
    """
    header = ['Family', 'Domain', 'Seed domains', 'Full region domains',
              'Rfam ID', 'Description', 'RNA type']

    # create a file for all families
    with open('domains/all-domains.csv', 'w') as f_out:
        csvwriter = csv.writer(f_out)
        csvwriter.writerow(header)
        for line in data:
            csvwriter.writerow(line)

    # Precompute parent group inclusion flags for each family
    parent_inclusion = {}
    for line in data:
        this_domain = line[1].lower()
        rfam_acc = line[0]
        full_domains = line[3]
        # Eukaryota inclusion logic
        is_eukaryote = (
            rfam_acc in WHITELIST or
            'eukaryota' in this_domain or
            ('fungi' in this_domain) or
            (get_fungi_percentage(full_domains) >= FUNGI_THRESHOLD)
        )
        parent_inclusion[rfam_acc] = {
            'Eukaryota': is_eukaryote,
            # Extend for other parent groups as needed
        }

    # create domain-specific files
    for domain in DOMAINS:
        if domain == 'Other':
            continue
        filename = 'domains/{}.csv'.format(domain.lower().replace(' ', '-'))
        with open(filename, 'w') as f_out:
            csvwriter = csv.writer(f_out)
            csvwriter.writerow(header)
            for line in data:
                this_domain = line[1].lower()
                rfam_acc = line[0]
                full_domains = line[3]
                if this_domain in ['bacteria/eukaryota']:
                    continue
                # Whitelist always included
                if rfam_acc in WHITELIST:
                    csvwriter.writerow(line)
                    continue
                # For subgroups: include if Domain field contains the subgroup name OR if threshold met with parent inclusion
                if domain in SUBGROUP_PARENT:
                    parent = SUBGROUP_PARENT[domain]
                    domain_field = line[1]
                    # Include if the domain name appears in the Domain field (handles Fungi, Eukaryota+Fungi, Mixed/Fungi, Fungi/Mixed, etc.)
                    if domain.lower() in domain_field.lower():
                        csvwriter.writerow(line)
                        continue
                    # Also include if threshold met and parent group inclusion is True
                    if domain == 'Fungi' and get_fungi_percentage(full_domains) >= FUNGI_THRESHOLD:
                        if parent_inclusion[rfam_acc].get(parent, False):
                            csvwriter.writerow(line)
                        continue
                # Standard inclusion for major domains
                if domain.lower() in this_domain:
                    csvwriter.writerow(line)
                    continue


def update_summary():
    """
    Update summary.md file with domain counts.
    """
    summary_file = 'domains/Readme.md'
    with open(summary_file, 'w') as f_out:
        header = """# Summary

The number of Rfam families observed in different domains:

```
"""
        f_out.write(header)
    cmd = ("cut -d ',' -f 2,2 domains/all-domains.csv | sort | uniq -c | "
           # "grep -v Mixed | "
           # "grep -v '/' | "
           "grep -v Domain | "
           "sort -nr >> {}".format(summary_file))
    os.system(cmd)
    os.system("echo '```' >> {}".format(summary_file))


@click.command()
@click.option('--precompute-seed', is_flag=True, required=False, help='Store seed data files')
@click.option('--precompute-full', is_flag=True, required=False, help='Store full data files')
@click.option('--cutoff', required=False, default=DOMAIN_CUTOFF, help='Percent of hits from the same domain')
def main(precompute_seed, precompute_full, cutoff):

    if precompute_seed:
        precompute_taxonomic_information('seed')
    if precompute_full:
        precompute_taxonomic_information('full')

    data = []
    for family in get_rfam_families():
        data.append(analyse_seed_full_taxonomic_distribution(family, cutoff))

    write_output_files(data)
    update_summary()


if __name__ == '__main__':
    main()
