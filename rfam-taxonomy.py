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

from scripts.rfam_db import get_rfam_families, get_taxonomy_info


DATA_SEED_PATH = 'data-seed'
DATA_FULL_REGION_PATH = 'data-full-region'
DOMAINS = sorted(['Archaea', 'Bacteria', 'Eukaryota', 'Other', 'unclassified sequences', 'Viruses', 'Viroids'])
DOMAIN_CUTOFF = 90 # at least 90% of sequences must be from this domain


def precompute_taxonomic_information(analysis_type):
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


def analyse_taxonomic_distribution(rfam_acc, DATA_PATH):
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
    for domain in DOMAINS:
        data[domain] = round(data[domain]*100.0/total, 2)
    return data


def get_major_domain(seed):
    major_domain = 'Mixed'
    maximum = max(seed, key=seed.get)
    if seed[maximum] >= DOMAIN_CUTOFF:
        major_domain = maximum
    return major_domain


def get_domains(data):
    output = []
    for domain, proportion in sorted(data.items(), key=lambda x: x[1], reverse=True):
        if proportion > 0:
            output.append('{} ({}%)'.format(domain, proportion))
    return ', '.join(output)



def analyse_seed_full_taxonomic_distribution(rfam_acc):
    seed = analyse_taxonomic_distribution(rfam_acc, DATA_SEED_PATH)
    full = analyse_taxonomic_distribution(rfam_acc, DATA_FULL_REGION_PATH)

    major_domain_seed = get_major_domain(seed)
    seed_domains = get_domains(seed)
    major_domain_full = get_major_domain(full)
    full_domains = get_domains(full)
    if major_domain_seed and major_domain_seed == major_domain_full:
        return [rfam_acc, major_domain_seed, seed_domains, full_domains]
    else:
        return [rfam_acc, '{}/{}'.format(major_domain_seed, major_domain_full), seed_domains, full_domains]


@click.command()
@click.option('--precompute-seed', is_flag=True, required=False, help='Store seed data files')
@click.option('--precompute-full', is_flag=True, required=False, help='Store full data files')
def main(precompute_seed, precompute_full):

    if precompute_seed:
        precompute_taxonomic_information('seed')
    if precompute_full:
        precompute_taxonomic_information('full')

    with open('domains.csv', 'w') as f_out:
        csvwriter = csv.writer(f_out)
        header = ['Family', 'Domain', 'Seed domains', 'Full region domains']
        csvwriter.writerow(header)
        for family in get_rfam_families():
            line = analyse_seed_full_taxonomic_distribution(family['rfam_acc'])
            csvwriter.writerow(line)

    cmd = ("cut -d ',' -f 2,2 domains.csv | sort | uniq -c | "
           "grep -v Mixed | "
           "grep -v '/' | "
           "grep -v Domain | "
           "sort -nr")
    os.system(cmd)


if __name__ == '__main__':
    main()
