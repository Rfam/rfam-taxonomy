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

from scripts.rfam_db import get_rfam_families, get_taxonomy_info, get_clan_membership


DATA_SEED_PATH = 'data-seed'
DATA_FULL_REGION_PATH = 'data-full-region'

# Domains (top-level biological domains)
DOMAINS = [
    'Archaea',
    'Bacteria',
    'Eukaryota',
    'Viruses',
    'Viroids',
    'unclassified sequences',
    'Other',
]

# Subgroups (taxonomic groups within parent domains)
SUBGROUPS = [
    # Archaea subgroups
    'Euryarchaeota',
    'Crenarchaeota',
    'Thaumarchaeota',
    # Bacteria subgroups
    'Proteobacteria',
    'Firmicutes',
    'Actinobacteria',
    'Bacteroidetes',
    'Cyanobacteria',
    # Eukaryota subgroups
    'Fungi',
    'Viridiplantae',
    'Metazoa',
    'Amoebozoa',
    'Alveolata',
    'Stramenopiles',
    # Virus subgroups
    'Flaviviridae',
    'Coronaviridae',
    'Retroviridae',
    'Herpesviridae',
    'Picornaviridae',
    'Orthomyxoviridae',
]

# All groups (combination of domains and subgroups)
ALL_GROUPS = DOMAINS + SUBGROUPS
DOMAIN_CUTOFF = 90 # at least 90% of sequences must be from this domain
# Ensure DOMAIN_CUTOFF is valid (51-100)
if not (51 <= DOMAIN_CUTOFF <= 100):
    raise ValueError(
        f"DOMAIN_CUTOFF ({DOMAIN_CUTOFF}) is invalid. It must be between 51 and 100 (inclusive) to ensure mutually exclusive major domains."
    )

SUBGROUP_THRESHOLD = 5  # include families with at least 5% of a subgroup in full regions
EPSILON = 1e-6 # small value used for floating-point comparison to avoid precision errors

# Mapping of subgroups to their parent domains
# To add a new subgroup, add it to the SUBGROUPS list above and to this mapping
SUBGROUP_PARENT = {
    # Archaea subgroups
    'Euryarchaeota': 'Archaea',
    'Crenarchaeota': 'Archaea',
    'Thaumarchaeota': 'Archaea',
    # Bacteria subgroups
    'Proteobacteria': 'Bacteria',
    'Firmicutes': 'Bacteria',
    'Actinobacteria': 'Bacteria',
    'Bacteroidetes': 'Bacteria',
    'Cyanobacteria': 'Bacteria',
    # Eukaryota subgroups
    'Fungi': 'Eukaryota',
    'Viridiplantae': 'Eukaryota',
    'Metazoa': 'Eukaryota',
    'Amoebozoa': 'Eukaryota',
    'Alveolata': 'Eukaryota',
    'Stramenopiles': 'Eukaryota',
    # Virus subgroups
    'Flaviviridae': 'Viruses',
    'Coronaviridae': 'Viruses',
    'Retroviridae': 'Viruses',
    'Herpesviridae': 'Viruses',
    'Picornaviridae': 'Viruses',
    'Orthomyxoviridae': 'Viruses',
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
    for domain in ALL_GROUPS:
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
                if taxon in ALL_GROUPS:
                    data[taxon] += count
                else:
                    data['Other'] += count
            total += count
    if num_rows == 0:
        return 'NO_DATA'
    if total != 0:
        for domain in ALL_GROUPS:
            data[domain] = round(data[domain]*100.0/total, 2)
    return data


def validate_multi_group_pattern(major_groups, cutoff):
    """
    Validate that A+B patterns only occur when groups have parent-child relationships.
    
    Raises ValueError if the pattern is invalid (e.g., multiple unrelated parent domains).
    """
    parents_in_set = [d for d in major_groups if d in SUBGROUP_PARENT.values()]
    subgroups_in_set = [d for d in major_groups if d in SUBGROUP_PARENT.keys()]
    
    # Valid A+B requires: exactly one parent, and all others are its subgroups
    if len(parents_in_set) == 1:
        parent = parents_in_set[0]
        for subgroup in subgroups_in_set:
            if SUBGROUP_PARENT[subgroup] != parent:
                raise ValueError(
                    f"Invalid A+B pattern: {'+'.join(sorted(major_groups))}. "
                    f"Subgroup '{subgroup}' is not a child of parent '{parent}'."
                )
    elif len(parents_in_set) > 1:
        raise ValueError(
            f"Invalid A+B pattern: {'+'.join(sorted(major_groups))}. "
            f"Multiple unrelated parent domains >= {cutoff}%: {parents_in_set}"
        )
    # If no parents (all subgroups), also invalid unless they share same parent
    elif len(parents_in_set) == 0 and len(subgroups_in_set) > 1:
        parent_set = set(SUBGROUP_PARENT[sg] for sg in subgroups_in_set)
        if len(parent_set) > 1:
            raise ValueError(
                f"Invalid A+B pattern: {'+'.join(sorted(major_groups))}. "
                f"Subgroups from different parents: {subgroups_in_set}"
            )


def get_major_group(data, cutoff):
    """
    Find the prevalent group (domain or subgroup, for example, Eukaryota or Fungi):

    {'Eukaryota': 100.0, 'Other': 0.0, 'Viruses': 0.0, 'unclassified sequences': 0.0, 'Viroids': 0.0, 'Archaea': 0.0, 'Bacteria': 0.0}
    """
    if data == 'NO_DATA':
        return 'No Data'
    
    # Remove unclassified sequences and renormalize if present
    if 'unclassified sequences' in data and data['unclassified sequences'] > 0:
        unclassified_pct = data['unclassified sequences']
        if unclassified_pct >= (100.0 - EPSILON):
            return 'Mixed'  # Only unclassified sequences present
        
        # Scale up remaining percentages to sum to 100%
        # Use multiplication by scale factor rather than summing percentages, because
        # parent+subgroup double-counting means percentages sum to >100% (e.g., a family
        # with 10% Fungi contributes to both Fungi and Eukaryota percentages).
        # Scaling preserves the proportional relationships while renormalizing.
        scale_factor = 100.0 / (100.0 - unclassified_pct)
        data = {
            domain: value * scale_factor
            for domain, value in data.items()
            if domain != 'unclassified sequences'
        }
    
    # Find all groups above cutoff
    major_groups = [domain for domain, value in data.items() if value >= cutoff]
    if len(major_groups) == 1:
        return major_groups[0]
    elif len(major_groups) > 1:
        validate_multi_group_pattern(major_groups, cutoff)
        return '+'.join(sorted(major_groups))
    else:
        return 'Mixed'


def get_groups(data):
    """
    List all groups (domains and subgroups) in which a family has been observed.
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



def get_group_percentage(groups_str, group, allow_prefix=True):
    """
    Extract the percentage for a given group (domain or subgroup) from a groups string.

    Args:
        groups_str: String with group percentages (e.g., 'Bacteria (0.25%), Eukaryota (40%)')
        group: Group name to search for
        allow_prefix: If True, matches "Parent:Group (X%)" format. If False, only standalone "Group (X%)"

    Examples:
        get_group_percentage('Eukaryota:Fungi (9.66%)', 'Fungi') => 9.66
        get_group_percentage('Eukaryota (40%)', 'Eukaryota') => 40.0
        get_group_percentage('Eukaryota:Fungi (9.66%)', 'Eukaryota', allow_prefix=False) => 0.0

    Returns:
        Percentage as float, or 0.0 if not found
    """
    if allow_prefix:
        # Match "group (x%)" possibly preceded by a parent prefix like "Eukaryota:Fungi"
        pattern = r'[^,]*\b{} \(([0-9.]+)%\)'.format(re.escape(group))
    else:
        # Only match standalone group, not as part of "Parent:Group" format
        pattern = r'\b{} \(([0-9.]+)%\)'.format(re.escape(group))
    
    match = re.search(pattern, groups_str)
    if match:
        return float(match.group(1))
    return 0.0


def analyse_seed_full_taxonomic_distribution(family, cutoff):
    """
    Compare groups observed in seed alignments and full region hits.
    """
    seed = get_taxonomic_distribution(family['rfam_acc'], DATA_SEED_PATH)
    full = get_taxonomic_distribution(family['rfam_acc'], DATA_FULL_REGION_PATH)

    major_group_seed = get_major_group(seed, cutoff)
    seed_groups = get_groups(seed)

    major_group_full = get_major_group(full, cutoff)
    full_groups = get_groups(full)

    # Validate: seed should never be 'No Data' (Rfam requires at least 2 sequences in seed)
    if major_group_seed == 'No Data':
        raise ValueError(
            f"Impossible state: seed has 'No Data' for {family['rfam_acc']}. "
            f"Rfam requires at least 2 sequences in seed alignments."
        )

    # Rule 1: If both are the same and not Mixed, use that value
    if major_group_seed == major_group_full and major_group_seed != 'Mixed':
        domain_field = major_group_seed
    # Rule 2: If either is Mixed or full is No Data, use the appropriate combination
    elif major_group_seed == 'Mixed' or major_group_full in ('Mixed', 'No Data'):
        # Handle all combinations (seed is never 'No Data' thanks to check above)
        if major_group_seed == 'Mixed' and major_group_full == 'Mixed':
            domain_field = 'Mixed'
        elif major_group_seed == 'Mixed' and major_group_full == 'No Data':
            domain_field = 'Mixed/No Data'
        elif major_group_full == 'No Data':
            domain_field = '{}/No Data'.format(major_group_seed)
        elif major_group_seed == 'Mixed':
            domain_field = 'Mixed/{}'.format(major_group_full)
        else:  # major_group_full == 'Mixed'
            domain_field = '{}/Mixed'.format(major_group_seed)
    else:
        # Rule 3: Handle standard A/B format (A or B could be C+D, so we can get A/C+D or C+D/B)
        domain_field = '{}/{}'.format(major_group_seed, major_group_full)
    return [
        family['rfam_acc'],
        domain_field,
        seed_groups,
        full_groups,
        family['rfam_id'],
        family['description'],
        family['type'],
    ]


def write_output_files(data):
    """
    Generate output files.
    """
    print("Generating domain-specific CSV files...")
    header = ['Family', 'Domain', 'Seed groups', 'Full region groups',
              'Rfam ID', 'Description', 'RNA type']

    # Write all families to all-domains.csv
    with open('domains/all-domains.csv', 'w') as f_out:
        csvwriter = csv.writer(f_out)
        csvwriter.writerow(header)
        for line in data:
            csvwriter.writerow(line)
    print("  Created domains/all-domains.csv")

    # Sort domains: subgroups first, then parents, then others
    # This ensures we know which families are in subgroups before writing parent files
    subgroups = SUBGROUPS
    # Use deterministic ordering for parent domains to keep CSV generation stable
    parents = sorted(set(SUBGROUP_PARENT.values()))
    other_domains = [d for d in DOMAINS if d not in parents and d != 'Other']
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
                full_groups = line[3]
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
                    seed_groups = line[2]
                    # (B1) Domain name appears in domain_field
                    if domain.lower() in domain_field.lower():
                        csvwriter.writerow(line)
                        subgroup_families[domain].add(rfam_acc)
                        continue
                    # (B2) Subgroup threshold: only if parent domain is substantively present in seed
                    # This ensures subgroups are restricted to the scope of their parent
                    # (prevents viral/bacterial contamination in eukaryotic subgroups like Fungi)
                    if get_group_percentage(full_groups, domain) >= SUBGROUP_THRESHOLD:
                        parent_pct_in_seed = get_group_percentage(seed_groups, parent, allow_prefix=False)
                        if parent_pct_in_seed >= SUBGROUP_THRESHOLD:
                            csvwriter.writerow(line)
                            subgroup_families[domain].add(rfam_acc)
                        continue
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
                
        print(f"  Created {filename}")
    
    print("Done generating CSV files")

def update_summary():
    """
    Update summary.md file with per-file statistics and group classification counts.
    """
    import glob
    summary_file = 'domains/Readme.md'
    
    # Count families and clans in each CSV/clanin file
    file_stats = {}
    for csv_file in glob.glob('domains/*.csv'):
        if csv_file.endswith('all-domains.csv'):
            continue
        
        # Count families (rows in CSV, excluding header)
        family_count = 0
        with open(csv_file, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # Skip header
            family_count = sum(1 for row in reader)
        
        # Count clans (lines in clanin file)
        clan_count = 0
        clanin_file = csv_file.replace('.csv', '.clanin')
        if os.path.exists(clanin_file):
            with open(clanin_file, 'r') as f:
                clan_count = sum(1 for line in f)
        
        # Extract domain/subgroup name from filename
        name = os.path.basename(csv_file).replace('.csv', '').replace('-', ' ')
        file_stats[name] = {'families': family_count, 'clans': clan_count}
    
    # Read all domain field values from all-domains.csv for pattern counts
    domain_counts = {}
    with open('domains/all-domains.csv', 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            if len(row) >= 2:
                domain_field = row[1]
                domain_counts[domain_field] = domain_counts.get(domain_field, 0) + 1
    
    # Write summary
    with open(summary_file, 'w') as f_out:
        f_out.write("""# Summary

```
The number of Rfam families and clans for each group (number of lines in .csv / .clanin files):
""")
        # Write file statistics sorted by family count (descending)
        sorted_files = sorted(file_stats.items(), key=lambda x: x[1]['families'], reverse=True)
        for name, stats in sorted_files:
            f_out.write(f"{stats['families']:7d} / {stats['clans']:<4d} {name}\n")
        
        f_out.write("\nDistribution of Domain field classifications:\n")
        f_out.write("(See ../Readme.md for explanation of Domain field format)\n")
        # Write all domain field values sorted by count
        sorted_domains = sorted(domain_counts.items(), key=lambda x: x[1], reverse=True)
        for domain_field, count in sorted_domains:
            f_out.write(f"{count:7d} {domain_field}\n")
        
        f_out.write("```\n")


def generate_clanin_files():
    """
    Generate domain-specific clanin files from Rfam database clan information.
    
    For each domain, creates a clanin file containing only clans with >1 family
    present in that domain. Queries the database directly to get current clan membership.
    """
    print("Retrieving clan membership from database...")
    clan_data = get_clan_membership()
    
    print("Generating domain-specific clanin files...")
    for domain in ALL_GROUPS:
        if domain == 'Other':
            continue
        
        csv_filename = 'domains/{}.csv'.format(domain.lower().replace(' ', '-'))
        clanin_filename = 'domains/{}.clanin'.format(domain.lower().replace(' ', '-'))
        
        if not os.path.exists(csv_filename):
            continue
        
        # Extract Rfam IDs from domain CSV (5th column)
        domain_rfam_ids = set()
        with open(csv_filename, newline='') as csvfile:
            reader = csv.reader(csvfile)
            next(reader)  # Skip header
            for row in reader:
                if len(row) >= 5:
                    domain_rfam_ids.add(row[4].strip())
                else:
                    print(f"Skipping malformed row in {csv_filename}: {row}")
        
        # Write clanin file for this domain
        with open(clanin_filename, 'w') as outfile:
            for clan_acc in clan_data.keys():
                clan_rfam_ids = clan_data[clan_acc]
                # Only keep Rfam IDs present in this domain
                present = [rfam_id for rfam_id in clan_rfam_ids if rfam_id in domain_rfam_ids]
                # Only write clans with >1 family in this domain
                if len(present) > 1:
                    outfile.write(f"{clan_acc} {' '.join(present)}\n")
        
        print(f"  Created {clanin_filename}")
    
    print("Done generating clanin files")


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
    
    # Generate domain-specific clanin files from database
    generate_clanin_files()


if __name__ == '__main__':
    main()
