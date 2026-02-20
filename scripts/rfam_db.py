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


import pymysql.cursors

# Connect to the public Rfam database
def get_connection():
    connection = pymysql.connect(host='mysql-rfam-public.ebi.ac.uk',
                                 user='rfamro',
                                 db='Rfam',
                                 port=4497,
                                 cursorclass=pymysql.cursors.DictCursor)
    return connection


def get_rfam_families():
    data = []
    try:
        connection = get_connection()
        with connection.cursor() as cursor:
            sql = """SELECT rfam_acc, rfam_id, description, type
                     FROM family
                     ORDER BY rfam_acc"""
            cursor.execute(sql)
            for result in cursor.fetchall():
                data.append(result)
    finally:
        connection.close()
    return data


def get_taxonomy_info(rfam_acc, analysis_type):

    if analysis_type == 'full':
        sql = """
        SELECT t3.tax_string, t3.ncbi_id, count(*) as count
        FROM full_region t1, rfamseq t2, taxonomy t3
        WHERE t1.rfamseq_acc = t2.rfamseq_acc
        AND is_significant = 1
        AND t2.ncbi_id = t3.ncbi_id
        AND rfam_acc = '{rfam_acc}'
        GROUP BY t3.tax_string
        """
    elif analysis_type == 'seed':
        sql = """
        SELECT t3.tax_string, t3.ncbi_id, count(*) as count
        FROM seed_region t1, rfamseq t2, taxonomy t3
        WHERE t1.rfamseq_acc = t2.rfamseq_acc
        AND t2.ncbi_id = t3.ncbi_id
        AND rfam_acc = '{rfam_acc}'
        GROUP BY t3.tax_string
        """
    data = []
    try:
        connection = get_connection()
        with connection.cursor() as cursor:
            cursor.execute(sql.format(rfam_acc=rfam_acc))
            for result in cursor.fetchall():
                data.append((result['tax_string'], result['count'], result['ncbi_id']))
    finally:
        connection.close()
    return data


def get_clan_membership():
    """
    Retrieve clan membership information from Rfam database.
    
    Returns:
        Dictionary mapping clan_acc to list of rfam_ids in that clan.
        Example: {'CL00001': ['5S_rRNA', 'tRNA', ...], 'CL00002': [...]}
    """
    clan_data = {}
    try:
        connection = get_connection()
        with connection.cursor() as cursor:
            sql = """
            SELECT c.clan_acc, f.rfam_id
            FROM clan c
            JOIN clan_membership cm ON c.clan_acc = cm.clan_acc
            JOIN family f ON cm.rfam_acc = f.rfam_acc
            ORDER BY c.clan_acc, f.rfam_id
            """
            cursor.execute(sql)
            for result in cursor.fetchall():
                clan_acc = result['clan_acc']
                rfam_id = result['rfam_id']
                if clan_acc not in clan_data:
                    clan_data[clan_acc] = []
                clan_data[clan_acc].append(rfam_id)
    finally:
        connection.close()
    return clan_data
