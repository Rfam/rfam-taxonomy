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
import click

from scripts.rfam_db import get_rfam_families, get_taxonomy_info, get_clan_membership


DATA_SEED_PATH = 'data-seed'
DATA_FULL_REGION_PATH = 'data-full-region'

DOMAINS = sorted([
    'Archaea',
    'Bacteria',
    'Eukaryota',
    'unclassified sequences',
    'Viruses',
    'Viroids',
    'Other',
])
DOMAIN_CUTOFF = 90 # at least 90% of sequences must be from this domain

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
    with open(os.path.join(DATA_PATH, '{}.csv'.format(rfam_acc)), 'r') as f:
        csvreader = csv.reader(f)
        for row in csvreader:
            tax_string, count, _ = row
            count = int(count)
            taxon = tax_string.split(';')[0]

            if '.' in taxon:
                taxon = taxon.replace('.', '') # example: Bacteria.
            if taxon == 'Unclassified':
                taxon = 'unclassified sequences'

            if taxon in DOMAINS:
                data[taxon] += count
            else:
                data['Other'] += count
            total += count
    if total != 0:
        for domain in DOMAINS:
            data[domain] = round(data[domain]*100.0/total, 2)
    return data


def get_major_domain(data, cutoff):
    """
    Find the prevalent domain (for example, Eukaryota):

    {'Eukaryota': 100.0, 'Other': 0.0, 'Viruses': 0.0, 'unclassified sequences': 0.0, 'Viroids': 0.0, 'Archaea': 0.0, 'Bacteria': 0.0}
    """
    major_domain = 'Mixed'
    maximum = max(data, key=data.get)
    if data[maximum] >= cutoff:
        major_domain = maximum
    else:
        # get distinct domains
        found_domains = []
        for domain, value in iter(data.items()):
            if value > 0:
                found_domains.append(domain)
        # if only two domains and one of them is `unclassified`, consider the other one major domain
        if len(found_domains) == 2 and 'unclassified sequences' in found_domains:
            found_domains.remove('unclassified sequences')
            major_domain = found_domains.pop()
    return major_domain


def get_domains(data):
    """
    List all domains in which a family has been observed.
    """
    output = []
    for domain, proportion in sorted(data.items(), key=lambda x: x[1], reverse=True):
        if proportion > 0:
            output.append('{} ({}%)'.format(domain, proportion))
    return ', '.join(output)


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

    if major_domain_seed and major_domain_seed == major_domain_full:
        return [
            family['rfam_acc'],
            major_domain_seed,
            seed_domains,
            full_domains,
            family['rfam_id'],
            family['description'],
            family['type'],
        ]
    else:
        return [
            family['rfam_acc'],
            '{}/{}'.format(major_domain_seed, major_domain_full),
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
    print("Generating domain-specific CSV files...")
    header = ['Family', 'Domain', 'Seed domains', 'Full region domains',
              'Rfam ID', 'Description', 'RNA type']

    # create a file for all families
    with open('domains/all-domains.csv', 'w') as f_out:
        csvwriter = csv.writer(f_out)
        csvwriter.writerow(header)
        for line in data:
            csvwriter.writerow(line)
    print("  Created domains/all-domains.csv")

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
                if this_domain in ['Bacteria/Eukaryota']: # skip families that have Bacteria in SEED but mostly Eukaryotes in full
                    continue
                elif rfam_acc in WHITELIST:
                    csvwriter.writerow(line)
                elif domain.lower() in this_domain:
                    csvwriter.writerow(line)
        print(f"  Created {filename}")
    print("Done generating CSV files")


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


def generate_clanin_files():
    """
    Generate domain-specific clanin files from Rfam database clan information.
    
    For each domain, creates a clanin file containing only clans with >1 family
    present in that domain. Queries the database directly to get current clan membership.
    """
    print("Retrieving clan membership from database...")
    clan_data = get_clan_membership()
    
    print("Generating domain-specific clanin files...")
    for domain in DOMAINS:
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
