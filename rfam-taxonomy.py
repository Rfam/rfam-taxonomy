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
SUBGROUP_THRESHOLD = 5  # include families with at least 5% of a subgroup in full regions

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

            # Check if any subgroup appears in the taxonomy lineage
            # Count in both the subgroup AND parent domain since subgroups are part of their parent
            for subgroup, parent in SUBGROUP_PARENT.items():
                if subgroup in tax_string:
                    data[subgroup] += count
                    data[parent] += count
                    break
            else:
                # No subgroup found, count in the top-level taxon
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


def validate_multi_domain_pattern(major_domains, cutoff):
    """
    Validate that A+B patterns only occur when domains have parent-child relationships.
    
    Raises ValueError if the pattern is invalid (e.g., multiple unrelated parent domains).
    """
    parents_in_set = [d for d in major_domains if d in SUBGROUP_PARENT.values()]
    subgroups_in_set = [d for d in major_domains if d in SUBGROUP_PARENT.keys()]
    
    # Valid A+B requires: exactly one parent, and all others are its subgroups
    if len(parents_in_set) == 1:
        parent = parents_in_set[0]
        for subgroup in subgroups_in_set:
            if SUBGROUP_PARENT[subgroup] != parent:
                raise ValueError(
                    f"Invalid A+B pattern: {'+'.join(sorted(major_domains))}. "
                    f"Subgroup '{subgroup}' is not a child of parent '{parent}'."
                )
    elif len(parents_in_set) > 1:
        raise ValueError(
            f"Invalid A+B pattern: {'+'.join(sorted(major_domains))}. "
            f"Multiple unrelated parent domains >= {cutoff}%: {parents_in_set}"
        )
    # If no parents (all subgroups), also invalid unless they share same parent
    elif len(parents_in_set) == 0 and len(subgroups_in_set) > 1:
        parent_set = set(SUBGROUP_PARENT[sg] for sg in subgroups_in_set)
        if len(parent_set) > 1:
            raise ValueError(
                f"Invalid A+B pattern: {'+'.join(sorted(major_domains))}. "
                f"Subgroups from different parents: {subgroups_in_set}"
            )


def get_major_domain(data, cutoff):
    """
    Find the prevalent domain (for example, Eukaryota):

    {'Eukaryota': 100.0, 'Other': 0.0, 'Viruses': 0.0, 'unclassified sequences': 0.0, 'Viroids': 0.0, 'Archaea': 0.0, 'Bacteria': 0.0}
    """
    if data == 'NO_DATA':
        return 'No Data'
    
    # Remove unclassified sequences and renormalize if present
    if 'unclassified sequences' in data and data['unclassified sequences'] > 0:
        # Calculate total percentage excluding unclassified
        total_excluding_unclassified = sum(
            value for domain, value in data.items() 
            if domain != 'unclassified sequences'
        )
        if total_excluding_unclassified == 0:
            return 'Mixed'  # Only unclassified sequences present
        # Renormalize to 100%
        data = {
            domain: (value / total_excluding_unclassified) * 100.0
            for domain, value in data.items()
            if domain != 'unclassified sequences'
        }
    
    # Find all domains above cutoff
    major_domains = [domain for domain, value in data.items() if value >= cutoff]
    if len(major_domains) == 1:
        return major_domains[0]
    elif len(major_domains) > 1:
        validate_multi_domain_pattern(major_domains, cutoff)
        return '+'.join(sorted(major_domains))
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

    # Subgroup percentages are already included in parent percentages
    # (computed in get_taxonomic_distribution)

    output = []
    for domain, proportion in sorted(data.items(), key=lambda x: x[1], reverse=True):
        if proportion > 0:
            # Show subgroups with their parent prefix
            display_name = domain
            if domain in SUBGROUP_PARENT:
                display_name = '{}:{}'.format(SUBGROUP_PARENT[domain], domain)
            output.append('{} ({:.2f}%)'.format(display_name, proportion))
    return ', '.join(output)



def get_domain_percentage(domains_str, domain, allow_prefix=True):
    """
    Extract the percentage for a given domain from a domains string.

    Args:
        domains_str: String with domain percentages (e.g., 'Bacteria (0.25%), Eukaryota (40%)')
        domain: Domain name to search for
        allow_prefix: If True, matches "Parent:Domain (X%)" format. If False, only standalone "Domain (X%)"

    Examples:
        get_domain_percentage('Eukaryota:Fungi (9.66%)', 'Fungi') => 9.66
        get_domain_percentage('Eukaryota (40%)', 'Eukaryota') => 40.0
        get_domain_percentage('Eukaryota:Fungi (9.66%)', 'Eukaryota', allow_prefix=False) => 0.0

    Returns:
        Percentage as float, or 0.0 if not found
    """
    if allow_prefix:
        # Match "domain (x%)" possibly preceded by a parent prefix like "Eukaryota:Fungi"
        pattern = r'[^,]*\b{} \(([0-9.]+)%\)'.format(re.escape(domain))
    else:
        # Only match standalone domain, not as part of "Parent:Domain" format
        pattern = r'\b{} \(([0-9.]+)%\)'.format(re.escape(domain))
    
    match = re.search(pattern, domains_str)
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

    # Validate: seed should never be 'No Data' (Rfam requires at least 2 sequences in seed)
    if major_domain_seed == 'No Data':
        raise ValueError(
            f"Impossible state: seed has 'No Data' for {family['rfam_acc']}. "
            f"Rfam requires at least 2 sequences in seed alignments."
        )

    # Rule 1: If both are the same and not Mixed, use that value
    if major_domain_seed == major_domain_full and major_domain_seed != 'Mixed':
        domain_field = major_domain_seed
    # Rule 2: If either is Mixed or full is No Data, use the appropriate combination
    elif major_domain_seed == 'Mixed' or major_domain_full in ('Mixed', 'No Data'):
        # Handle all combinations (seed is never 'No Data' thanks to check above)
        if major_domain_seed == 'Mixed' and major_domain_full == 'Mixed':
            domain_field = 'Mixed'
        elif major_domain_seed == 'Mixed' and major_domain_full == 'No Data':
            domain_field = 'Mixed/No Data'
        elif major_domain_full == 'No Data':
            domain_field = '{}/No Data'.format(major_domain_seed)
        elif major_domain_seed == 'Mixed':
            domain_field = 'Mixed/{}'.format(major_domain_full)
        else:  # major_domain_full == 'Mixed'
            domain_field = '{}/Mixed'.format(major_domain_seed)
    else:
        # Rule 3: Handle standard A/B format (A or B could be C+D, so we can get A/C+D or C+D/B)
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

    # Write all families to all-domains.csv
    with open('domains/all-domains.csv', 'w') as f_out:
        csvwriter = csv.writer(f_out)
        csvwriter.writerow(header)
        for line in data:
            csvwriter.writerow(line)

    # Sort domains: subgroups first, then parents, then others
    # This ensures we know which families are in subgroups before writing parent files
    subgroups = list(SUBGROUP_PARENT.keys())
    parents = list(set(SUBGROUP_PARENT.values()))
    other_domains = [d for d in DOMAINS if d not in subgroups and d not in parents and d != 'Other']
    sorted_domains = subgroups + parents + other_domains
    
    # Track which families are written to each subgroup file
    subgroup_families = {subgroup: set() for subgroup in subgroups}

    # Write domain-specific files
    for domain in sorted_domains:
        filename = 'domains/{}.csv'.format(domain.lower().replace(' ', '-'))
        with open(filename, 'w') as f_out:
            csvwriter = csv.writer(f_out)
            csvwriter.writerow(header)
            for line in data:
                this_domain = line[1].lower()
                rfam_acc = line[0]
                full_domains = line[3]
                domain_field = line[1]
                # Exclude bacteria/eukaryota special case
                if this_domain == 'bacteria/eukaryota':
                    continue
                # (A) Whitelist: always include
                if rfam_acc in WHITELIST:
                    csvwriter.writerow(line)
                    if domain in subgroups:
                        subgroup_families[domain].add(rfam_acc)
                    continue
                # (B) Subgroup file logic
                if domain in SUBGROUP_PARENT:
                    parent = SUBGROUP_PARENT[domain]
                    seed_domains = line[2]
                    # (B1) Domain name appears in domain_field
                    if domain.lower() in domain_field.lower():
                        csvwriter.writerow(line)
                        subgroup_families[domain].add(rfam_acc)
                        continue
                    # (B2) Subgroup threshold: only if parent domain is substantively present in seed
                    # This ensures subgroups are restricted to the scope of their parent
                    # (prevents viral/bacterial contamination in eukaryotic subgroups like Fungi)
                    if get_domain_percentage(full_domains, domain) >= SUBGROUP_THRESHOLD:
                        parent_pct_in_seed = get_domain_percentage(seed_domains, parent, allow_prefix=False)
                        if parent_pct_in_seed >= SUBGROUP_THRESHOLD:
                            csvwriter.writerow(line)
                            subgroup_families[domain].add(rfam_acc)
                        continue
                # (C) Parent domain file logic: include all families from child subgroups
                if domain in parents:
                    # Include families from any subgroup that has this domain as parent
                    for subgroup, parent in SUBGROUP_PARENT.items():
                        if parent == domain and rfam_acc in subgroup_families[subgroup]:
                            csvwriter.writerow(line)
                            break
                    else:
                        # Also check standard inclusion if not already written
                        if domain.lower() in this_domain:
                            csvwriter.writerow(line)
                    continue
                # (D) Standard inclusion: domain name in this_domain
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
