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

SUBGROUP_THRESHOLD = 5  # include families with at least 5% of a group in full regions (and its parent in seed)
# Also include families found in at least PRESENCE_THRESHOLD% of the taxa of a group that have any
# full region hit, and in at least PRESENCE_MIN_TAXA of them (and its parent in seed). This catches
# RNAs present throughout a group whose hits are outnumbered by other groups (e.g. U6 snRNA in Fungi:
# 66% of fungal taxa, but under 2% of the hits, because animal and plant genomes carry many copies).
PRESENCE_THRESHOLD = 5
PRESENCE_MIN_TAXA = 3
# Also include families with at least SEED_THRESHOLD% of their seed sequences from a group, and at
# least one; 0 means any one seed sequence. This catches RNAs restricted to a few lineages of a group
# whose seed includes them (e.g. the yeast snoRNAs snoZ122 and snoZ196: about a quarter of the seed is
# fungal, but under 1.5% of the full region hits).
SEED_THRESHOLD = 0
# The rules above apply only if at least PARENT_SEED_THRESHOLD% of the seed sequences, and at least
# one, are from the group's parent domain (for a domain, the domain itself); 0 means any one seed
# sequence. This keeps out families whose hits in a group are probably contamination, e.g. a
# bacterial or viral family with full region hits in eukaryotic genome assemblies.
PARENT_SEED_THRESHOLD = 0

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

# Current NCBI names (and older spellings) for the subgroups above, mapped to the name used in SUBGROUPS.
# NCBI Taxonomy renamed prokaryotic phyla in 2021-2024; Rfam lineages contain both forms.
SUBGROUP_SYNONYMS = {
    'Pseudomonadota': 'Proteobacteria',
    'Bacillota': 'Firmicutes',
    'Actinomycetota': 'Actinobacteria',
    'Bacteroidota': 'Bacteroidetes',
    'Cyanobacteriota': 'Cyanobacteria',
    'Thermoproteota': 'Crenarchaeota',
    'Nitrososphaerota': 'Thaumarchaeota',
    'Methanobacteriota': 'Euryarchaeota',
    'Halobacteriota': 'Euryarchaeota',
    'Thermoplasmatota': 'Euryarchaeota',
    'Candidatus Thermoplasmatota': 'Euryarchaeota',
    'stramenopiles': 'Stramenopiles',
    'Orthoherpesviridae': 'Herpesviridae',
}

WHITELIST = [
    'RF00001', # 5S
    'RF00005', # tRNA
    'RF01852', # tRNA-Sec
]

# Hand-curated additions to and removals from specific group files, for families the rules above
# get wrong. Both lists are checked against the data on every run (see print_curation_report()):
# re-check them at each Rfam release.
# Stand-ins for organellar rRNA, added to every eukaryotic subgroup below. Rfam has no mitochondrial
# or plastid rRNA families, and its bacterial SSU and LSU rRNA models are the ones that annotate
# organellar rRNA (in Rfam 15.1 they hit about 90% of the fungal and animal mitochondrial genomes with
# any Rfam hit, and most plastid genomes). Remove them once Rfam has organellar rRNA families.
ORGANELLAR_RRNA_STAND_INS = [
    'RF00177', # SSU_rRNA_bacteria: stand-in for mitochondrial and plastid SSU rRNA
    'RF02541', # LSU_rRNA_bacteria: stand-in for mitochondrial and plastid LSU rRNA
]
# Families always included in a group file (and so in its parent's file)
GROUP_INCLUDE = {
    'Fungi': [
        'RF00049', # SNORD36: strong hits in Saccharomycotina, Zoopagomycota and Mucoromycota
        'RF00563', # SNORA53: hits in 18 fungal lineages, mostly Mucoromycota; seed is all animal
    ] + ORGANELLAR_RRNA_STAND_INS,
    'Viridiplantae': ORGANELLAR_RRNA_STAND_INS,
    'Metazoa': ORGANELLAR_RRNA_STAND_INS,
    'Amoebozoa': ORGANELLAR_RRNA_STAND_INS,
    'Alveolata': ORGANELLAR_RRNA_STAND_INS,
    'Stramenopiles': ORGANELLAR_RRNA_STAND_INS,
}
# Families never included in a group file (checked last, overrides every other rule)
GROUP_EXCLUDE = {
    'Viridiplantae': [
        'RF00254', # mir-16: animal miRNA; plant hits are weak (best E-value 9e-5), near the GA threshold
        'RF02253', # IRE_II: animal iron-response element; plant hits are weak (best 36.4 bits, GA 30.0)
    ],
    'Alveolata': [
        'RF00254', # mir-16: animal miRNA; alveolate hits are weak (best E-value 0.02), near the GA threshold
    ],
}


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


def get_lineage_groups(tax_string):
    """
    Return the top-level taxon of a lineage and the groups its hits are counted in:
    a subgroup and its parent domain if a subgroup appears in the lineage, otherwise
    the top-level taxon (or 'Other' if it is not one of the groups).

    Example:
    'Eukaryota; Fungi; Dikarya; Ascomycota; ...' => ('Eukaryota', ['Fungi', 'Eukaryota'])
    """
    taxon = tax_string.split(';')[0]
    if '.' in taxon:
        taxon = taxon.replace('.', '') # example: Bacteria.
    if taxon == 'Unclassified':
        taxon = 'unclassified sequences'
    if taxon not in DOMAINS:  # a few lineages are stored from species to domain
        taxon = tax_string.rstrip('. ').split(';')[-1].strip()

    # Check if any subgroup appears in the taxonomy lineage
    # Count in both the subgroup AND parent domain since subgroups are part of their parent
    # Parse lineage into tokens to avoid false positives from substring matching;
    # strip trailing periods (e.g. "Bacteria; Firmicutes.") and map current NCBI names to subgroup names
    tax_tokens = {SUBGROUP_SYNONYMS.get(t, t)
                  for t in (token.strip().rstrip('.') for token in tax_string.split(';'))}
    for subgroup, parent in SUBGROUP_PARENT.items():
        if subgroup in tax_tokens:
            return taxon, [subgroup, parent]
    # No subgroup found, count in the top-level taxon
    if taxon in ALL_GROUPS:
        return taxon, [taxon]
    return taxon, ['Other']


def get_taxonomic_distribution(rfam_acc, DATA_PATH, exclude_unclassified=False):
    """
    Calculate the percentage of hits from each group (domain and subgroup) for a family.
    Subgroup hits are also counted in their parent domain, so percentages can sum to more than 100.

    If exclude_unclassified is True, rows from unclassified sequences are skipped.

    Example (groups with 0.0 omitted):
    {'Eukaryota': 45.51, 'Fungi': 9.66, 'Metazoa': 35.85, 'Bacteria': 48.6, 'Proteobacteria': 40.2, 'Archaea': 5.9}
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
            taxon, groups = get_lineage_groups(tax_string)
            if exclude_unclassified and taxon == 'unclassified sequences':
                continue
            for group in groups:
                data[group] += count
            total += count
    if num_rows == 0 or total == 0:
        return 'NO_DATA'
    if total != 0:
        for domain in ALL_GROUPS:
            data[domain] = round(data[domain]*100.0/total, 2)
    return data


def get_group_taxa(rfam_acc):
    """
    Return the set of taxa with full region hits in each group (domain and subgroup) for a family.
    Each row of the data files is one lineage string, which stops at the genus in most cases, and its
    NCBI taxid is one representative of that lineage (the query groups hits by lineage string). So the
    taxa counted here are lineages, roughly genera, identified by that representative taxid.

    Example:
    {'Eukaryota': {'4932', '9606', ...}, 'Fungi': {'4932', ...}, 'Metazoa': {'9606', ...}}
    """
    taxa = {}
    with open(os.path.join(DATA_FULL_REGION_PATH, '{}.csv'.format(rfam_acc)), 'r') as f:
        for tax_string, _, taxid in csv.reader(f):
            for group in get_lineage_groups(tax_string)[1]:
                taxa.setdefault(group, set()).add(taxid)
    return taxa


def get_group_seed_counts(rfam_acc):
    """
    Return the number of seed sequences from each group (domain and subgroup) for a family, and the
    total number of seed sequences. Subgroup sequences are also counted in their parent domain.

    Example:
    ({'Eukaryota': 23, 'Fungi': 5, 'Metazoa': 18}, 23)
    """
    counts = {}
    total = 0
    with open(os.path.join(DATA_SEED_PATH, '{}.csv'.format(rfam_acc)), 'r') as f:
        for tax_string, count, _ in csv.reader(f):
            for group in get_lineage_groups(tax_string)[1]:
                counts[group] = counts.get(group, 0) + int(count)
            total += int(count)
    return counts, total


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
    Find the prevalent group(s): every group with at least `cutoff` percent of the hits.
    Returns e.g. 'Eukaryota', 'Eukaryota+Fungi' (parent and subgroup both >= cutoff), 'Mixed', or 'No Data'.

    Example: {'Eukaryota': 100.0, 'Fungi': 95.0, 'Metazoa': 5.0, ...} => 'Eukaryota+Fungi'
    """
    if data == 'NO_DATA':
        return 'No Data'
    
    # Unclassified sequences carry no taxonomic signal: if they reach the cutoff themselves,
    # that is the answer; otherwise remove them and renormalize before looking for major groups
    unclassified_pct = data.get('unclassified sequences', 0)
    if unclassified_pct >= cutoff:
        return 'unclassified sequences'
    if unclassified_pct > 0:
        # Use multiplication by scale factor rather than summing percentages, because
        # parent+subgroup double-counting means percentages sum to >100% (e.g., a family
        # with 10% Fungi contributes to both Fungi and Eukaryota percentages).
        scale_factor = 100.0 / (100.0 - unclassified_pct)
        data = {
            domain: value * scale_factor
            for domain, value in data.items()
            if domain != 'unclassified sequences'
        }

    major_groups = [domain for domain, value in data.items() if value >= cutoff]
    if len(major_groups) == 1:
        return major_groups[0]
    elif len(major_groups) > 1:
        validate_multi_group_pattern(major_groups, cutoff)
        # Order with parent domains before subgroups to maintain Parent+Subgroup convention
        parents = [g for g in major_groups if g not in SUBGROUP_PARENT]
        subgroups = [g for g in major_groups if g in SUBGROUP_PARENT]
        return '+'.join(sorted(parents) + sorted(subgroups))

    # No group meets the cutoff
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



def get_group_percentage(groups_str, group):
    """
    Extract the percentage for a given group (domain or subgroup) from a groups string.

    Args:
        groups_str: String with group percentages (e.g., 'Eukaryota (40.00%), Eukaryota:Fungi (9.66%)')
        group: Group name to search for

    Examples:
        get_group_percentage('Eukaryota (40.00%), Eukaryota:Fungi (9.66%)', 'Fungi') => 9.66
        get_group_percentage('Eukaryota (40.00%), Eukaryota:Fungi (9.66%)', 'Eukaryota') => 40.0

    A parent name never matches the "Parent:" prefix of a subgroup entry, because the
    pattern requires the name to be followed by " (".

    Returns:
        Percentage as float, or 0.0 if not found
    """
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


def in_labels(group, labels):
    """
    True if the group name appears in any of the placement labels (e.g. 'Bacteria' in 'Bacteria+Firmicutes/Mixed').
    """
    return any(group.lower() in label.lower() for label in labels)


def get_placement_labels(line, cutoff):
    """
    Return the labels used to decide which group files a family goes into.

    Usually this is just the Domain field. Two exceptions:
    - Bacteria in seed but Eukaryota in full: the full hits are bacterial RNAs in eukaryotic
      assemblies (e.g. Wolbachia sRNAs in arthropod genomes), so only the seed label is used.
    - Seed mostly unclassified sequences (e.g. metagenome-derived seeds): the seed sequences
      with a known organism are classified with the same cutoff, and that label is added.
    """
    rfam_acc, domain_field = line[0], line[1]
    seed_part, _, full_part = domain_field.partition('/')
    if seed_part.split('+')[0] == 'Bacteria' and full_part.split('+')[0] == 'Eukaryota':
        return [seed_part]
    labels = [domain_field]
    if seed_part == 'unclassified sequences':
        classified = get_taxonomic_distribution(rfam_acc, DATA_SEED_PATH, exclude_unclassified=True)
        classified_group = get_major_group(classified, cutoff)
        if classified_group not in ('No Data', 'Mixed'):
            labels.append(classified_group)
    return labels


def write_output_files(data, cutoff, subgroup_threshold=SUBGROUP_THRESHOLD,
                       presence_threshold=PRESENCE_THRESHOLD, presence_min_taxa=PRESENCE_MIN_TAXA,
                       seed_threshold=SEED_THRESHOLD, parent_seed_threshold=PARENT_SEED_THRESHOLD):
    """
    Generate output files. The thresholds are explained where their defaults are defined.
    """
    print("Generating group-specific CSV files...")
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
    placement_labels = {line[0]: get_placement_labels(line, cutoff) for line in data}

    # Taxa with full region hits for each family and group, and in total for each group
    family_taxa = {line[0]: get_group_taxa(line[0]) for line in data}
    group_taxa_count = {}
    for group in ALL_GROUPS:
        group_taxa_count[group] = len(set().union(*(taxa.get(group, set()) for taxa in family_taxa.values())))

    def present_in(rfam_acc, group):
        """True if the family has full region hits in enough of the group's taxa (see PRESENCE_THRESHOLD)."""
        n = len(family_taxa[rfam_acc].get(group, ()))
        return n >= presence_min_taxa and n * 100.0 / group_taxa_count[group] >= presence_threshold

    # Seed sequences from each group, for each family
    family_seed = {line[0]: get_group_seed_counts(line[0]) for line in data}

    def in_seed(rfam_acc, group, threshold):
        """True if at least one seed sequence, and at least threshold percent of them, are from the group."""
        counts, total = family_seed[rfam_acc]
        n = counts.get(group, 0)
        return n >= 1 and n * 100.0 / total >= threshold

    def data_rules_include(line, domain, child_families):
        """
        Rules (B)-(D): True if the data place the family in the group file.
        For a parent domain, child_families gives the families in each subgroup file.
        """
        rfam_acc = line[0]
        full_groups = line[3]
        labels = placement_labels[rfam_acc]
        # (B) Subgroup file logic
        if domain in SUBGROUP_PARENT:
            parent = SUBGROUP_PARENT[domain]
            # (B1) Subgroup name appears in a placement label, or
            # (B2) subgroup threshold (share of hits), presence threshold (share of taxa) or
            # seed threshold (share of seed sequences): only if the parent domain is present in seed
            # This ensures subgroups are restricted to the scope of their parent
            # (prevents viral/bacterial contamination in eukaryotic subgroups like Fungi)
            return (in_labels(domain, labels) or
                    ((get_group_percentage(full_groups, domain) >= subgroup_threshold or
                      present_in(rfam_acc, domain) or
                      in_seed(rfam_acc, domain, seed_threshold)) and
                     in_seed(rfam_acc, parent, parent_seed_threshold)))
        # (C) Parent domain file logic: the same rules as for subgroups, plus every
        # family in a child subgroup file (so each subgroup file is a subset of its parent's)
        if domain in parents:
            return (in_labels(domain, labels) or
                    ((get_group_percentage(full_groups, domain) >= subgroup_threshold or
                      present_in(rfam_acc, domain) or
                      in_seed(rfam_acc, domain, seed_threshold)) and
                     in_seed(rfam_acc, domain, parent_seed_threshold)) or
                    any(rfam_acc in child_families[subgroup]
                        for subgroup, parent in SUBGROUP_PARENT.items() if parent == domain))
        # (D) Standard inclusion: domain name in a placement label
        return in_labels(domain, labels)

    # Track which families are written to each group file, and which ones the
    # data rules alone (without WHITELIST, GROUP_INCLUDE and GROUP_EXCLUDE) would include
    group_families = {domain: set() for domain in sorted_domains}
    data_families = {domain: set() for domain in sorted_domains}

    # Write domain-specific files
    for domain in sorted_domains:
        filename = 'domains/{}.csv'.format(domain.lower().replace(' ', '-'))
        with open(filename, 'w') as f_out:
            csvwriter = csv.writer(f_out)
            csvwriter.writerow(header)
            for line in data:
                rfam_acc = line[0]
                if data_rules_include(line, domain, data_families):
                    data_families[domain].add(rfam_acc)
                # (A) Whitelist: include only in Archaea, Bacteria, Eukaryota and their subgroups
                if rfam_acc in WHITELIST:
                    include = SUBGROUP_PARENT.get(domain, domain) in ['Archaea', 'Bacteria', 'Eukaryota']
                # (A2) Hand-curated additions (a parent gets them through rule (C))
                elif rfam_acc in GROUP_INCLUDE.get(domain, ()):
                    include = True
                # (B)-(D) Data rules
                else:
                    include = data_rules_include(line, domain, group_families)
                # (E) Hand-curated removals, checked last
                if rfam_acc in GROUP_EXCLUDE.get(domain, ()):
                    include = False
                if include:
                    csvwriter.writerow(line)
                    group_families[domain].add(rfam_acc)

        print(f"  Created {filename}")

    # Every family in a subgroup file must also be in its parent's file
    for subgroup, parent in SUBGROUP_PARENT.items():
        missing = group_families[subgroup] - group_families[parent]
        if missing:
            raise ValueError(f"Families in {subgroup} but not in its parent {parent}: {sorted(missing)}")

    print("Done generating CSV files")
    print_curation_report(data_families, family_taxa, group_taxa_count,
                          {line[0]: line[4] for line in data})


def check_curated_lists(rfam_accs):
    """
    Check GROUP_INCLUDE and GROUP_EXCLUDE: groups must exist, accessions must be Rfam families,
    and a family cannot be both included in and excluded from the same group.
    """
    for name, curated in (('GROUP_INCLUDE', GROUP_INCLUDE), ('GROUP_EXCLUDE', GROUP_EXCLUDE)):
        for group, accs in curated.items():
            if group not in ALL_GROUPS:
                raise ValueError(f"{name}: unknown group '{group}'")
            unknown = sorted(set(accs) - set(rfam_accs))
            if unknown:
                raise ValueError(f"{name}['{group}']: not Rfam families: {unknown}")
    for group in set(GROUP_INCLUDE) & set(GROUP_EXCLUDE):
        both = sorted(set(GROUP_INCLUDE[group]) & set(GROUP_EXCLUDE[group]))
        if both:
            raise ValueError(f"Families in both GROUP_INCLUDE and GROUP_EXCLUDE for {group}: {both}")


def print_curation_report(data_families, family_taxa, group_taxa_count, rfam_ids):
    """
    Report whether each GROUP_INCLUDE and GROUP_EXCLUDE entry is still needed: an entry is
    redundant if the data rules alone already make the same decision. The number and share of
    the group's taxa with full region hits are shown, so large changes between releases stand out.
    """
    print("Hand-curated group lists (GROUP_INCLUDE / GROUP_EXCLUDE):")
    for name, curated, needed_if in (('include', GROUP_INCLUDE, False), ('exclude', GROUP_EXCLUDE, True)):
        for group, accs in curated.items():
            for rfam_acc in accs:
                by_data = rfam_acc in data_families[group]
                n = len(family_taxa[rfam_acc].get(group, ()))
                status = 'needed' if by_data == needed_if else 'redundant, can be removed'
                print(f"  {name} {rfam_acc} ({rfam_ids[rfam_acc]}) {group}: hits in {n} of {group_taxa_count[group]} taxa "
                      f"({n * 100.0 / group_taxa_count[group]:.1f}%); {status}")

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
        sorted_files = sorted(file_stats.items(), key=lambda x: (-x[1]['families'], x[0]))
        for name, stats in sorted_files:
            f_out.write(f"{stats['families']:7d} / {stats['clans']:<4d} {name}\n")
        
        f_out.write("\nDistribution of Domain field classifications:\n")
        f_out.write("(See ../Readme.md for explanation of Domain field format)\n")
        # Write all domain field values sorted by count
        sorted_domains = sorted(domain_counts.items(), key=lambda x: (-x[1], x[0]))
        for domain_field, count in sorted_domains:
            f_out.write(f"{count:7d} {domain_field}\n")
        
        f_out.write("```\n")


def generate_clanin_files():
    """
    Generate group-specific clanin files from Rfam database clan information.
    
    For each group, creates a clanin file containing only clans with >1 family
    present in that group. Queries the database directly to get current clan membership.
    """
    print("Retrieving clan membership from database...")
    clan_data = get_clan_membership()
    
    print("Generating group-specific clanin files...")
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


@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--precompute-seed', is_flag=True, required=False, help='Store seed data files')
@click.option('--precompute-full', is_flag=True, required=False, help='Store full data files')
@click.option('--cutoff', required=False, default=DOMAIN_CUTOFF, type=click.IntRange(51, 100),
              help='Minimum percent of hits from one group (51-100)')
@click.option('--subgroup-threshold', required=False, default=SUBGROUP_THRESHOLD,
              type=click.FloatRange(0, 100),
              help='Include a family in a group file if at least this percent of its full region hits '
                   f'are from the group (0-100, default {SUBGROUP_THRESHOLD})')
@click.option('--presence-threshold', required=False, default=PRESENCE_THRESHOLD,
              type=click.FloatRange(0, 100),
              help='Include a family in a group file if it has full region hits in at least this percent '
                   f'of the group\'s taxa, and in at least --presence-min-taxa of them (0-100, '
                   f'default {PRESENCE_THRESHOLD})')
@click.option('--presence-min-taxa', required=False, default=PRESENCE_MIN_TAXA, type=click.IntRange(1),
              help=f'Minimum number of the group\'s taxa for --presence-threshold (default {PRESENCE_MIN_TAXA})')
@click.option('--seed-threshold', required=False, default=SEED_THRESHOLD, type=click.FloatRange(0, 100),
              help='Include a family in a group file if at least this percent of its seed sequences, '
                   f'and at least one, are from the group (0-100, default {SEED_THRESHOLD}: any one)')
@click.option('--parent-seed-threshold', required=False, default=PARENT_SEED_THRESHOLD,
              type=click.FloatRange(0, 100),
              help='Apply the three rules above (--subgroup-threshold, --presence-threshold, --seed-threshold) only if at least this percent of the seed sequences, and '
                   'at least one, are from the parent domain of the group (for a domain, the domain '
                   f'itself) (0-100, default {PARENT_SEED_THRESHOLD}: any one)')
def main(precompute_seed, precompute_full, cutoff, subgroup_threshold, presence_threshold,
         presence_min_taxa, seed_threshold, parent_seed_threshold):
    if precompute_seed:
        precompute_taxonomic_information('seed')
    if precompute_full:
        precompute_taxonomic_information('full')

    families = get_rfam_families()
    check_curated_lists([family['rfam_acc'] for family in families])

    data = []
    for family in families:
        data.append(analyse_seed_full_taxonomic_distribution(family, cutoff))

    write_output_files(data, cutoff, subgroup_threshold, presence_threshold, presence_min_taxa,
                       seed_threshold, parent_seed_threshold)

    # Generate domain-specific clanin files from database
    generate_clanin_files()

    # Summary counts lines in the .clanin files, so it must run after generate_clanin_files()
    update_summary()


if __name__ == '__main__':
    main()
