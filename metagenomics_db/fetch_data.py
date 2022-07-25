#-*- coding: utf-8 -*-

# Copyright 2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import gzip
import hashlib
import os
import sys
import re
import shutil
from os import stat

import pandas as pd
from Bio import SeqIO

def check_study_accession(study_accession):
    study_accession_re = re.compile(r"((E|D|S)RP[0-9]{6,})")
    if not study_accession or not study_accession_re.match(study_accession):
        logging.error(f"{study_accession} is not an valid Secondary study accession.")
        sys.exit(1)
    else:
        print("Starting the workflow")

def parsing_header(header, pipeline_version='5.0'):
    callers_versions = {
        '5.0': {'FGS': '1.31', 'Prodigal': '2.6.3'},
        '4.1': {'FGS': '1.20', 'Prodigal': '2.6.3'}}  # pipeline version and corresponding versions of callers
    partial_dict = {'00': 'full', '01': 'partial', '10': 'partial', '11': 'partial'}

#     length, coverage = get_length_and_coverage(header)
    if 'partial' in header:
        # prodigal header:
        # v5.0
        # >ERZ1738565.7-k99-1872_1 # 1 # 426 # -1 # ID=1_1;partial=10;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.615
        # v4.1
        # >ERZ476429.1-NODE-1-length-477407-cov-27.154_1 # 3 # 296 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.449 IPR017896/PF12838/8-78|G3DSA:3.30.70.20/18-69;2-17;70-90

        caller = 'Prodigal'
        start_coordinate, stop_coordinate, strand, *_ = re.findall(r"#\s(.*?)\s", header)
        id, partial, start_type, stop_type, rbs_motif, *_ = re.findall(r"\=(.*?)\;", header)

        # Flip strand if prediction is on the reverse strand
        if strand == '-1':
            partial = partial[::-1]
    else:
        # FGS header:
        # v5.0
        # >ERZ476429.1-NODE-1-length-477407-cov-27.154_120011_120250_+
        # v4.1
        # >ERZ476429.1-NODE-1-length-477407-cov-27.154_120011_120250_+ PF13412/19-65|IPR011991/G3DSA:1.10.10.10/15-78

        partial = '11'
        caller = 'FGS'
        list_fields = header.split(' ')[0].split('_')
        length = len(list_fields)
        sign_strand, stop_coordinate, start_coordinate = list_fields[length-1], list_fields[length-2], list_fields[length-3]
        strand = str(int(sign_strand + '1'))
    caller_field = caller + '_' + callers_versions[str(pipeline_version)][caller]
    partial_field = partial_dict[partial]
    return partial_field, partial, start_coordinate, stop_coordinate, strand, caller_field

def create_digest(input_string): # input_string=partial_+record.seq
    """
        Create digest sha128
        :param input_string: input string
        :return: digest
        """
    digest = hashlib.sha256(str(input_string).encode('utf-8')).hexdigest()
    return digest



def rename_contigs(assembly: str, out_folder:str,sequence_dir: str ):
    with gzip.open(os.path.join(out_folder, assembly+".faa.gz"), 'wt') as fasta_out,open(os.path.join(out_folder, assembly+"_contig_info.txt"), 'w') as peptides_metadata:
        with gzip.open(os.path.join(sequence_dir, assembly+"_FASTA_predicted_cds.faa.gz"), 'rt') as cds_fasta:
            for record in SeqIO.parse(cds_fasta, "fasta"):
                contig_name=record.id
                sequence = record.seq
                partial_field, partial, start_coordinate, stop_coordinate, strand, caller = parsing_header(record.description)
                digest = create_digest(contig_name + sequence)
                peptides_metadata.write('\t'.join([assembly, contig_name, digest, partial_field, start_coordinate, stop_coordinate, strand, caller]) + '\n')
                record.id=digest
                record.description=digest
                SeqIO.write(record, fasta_out, "fasta")




def count_proteins(sample_assembly_dict: dict, seq_dir: str):
    protein_counts = {}
    for sample, assemblies in sample_assembly_dict.items():
        total_count=0
        for assembly in assemblies:
            assembly_path=os.path.join(seq_dir,assembly+"_FASTA_predicted_cds.faa.gz")
            with gzip.open(assembly_path, 'rt') as assembly_fasta:
                for line in assembly_fasta:
                    if line.startswith(">"):
                        total_count += 1
        protein_counts[sample]=total_count
    return protein_counts

def build_db(study: str, db_dir:str,assembly_dir: str, cluster_dict: dict):
    for group_no, assembly_list in cluster_dict.items():
        total_count=0
        db_name = study + "_cluster_set_" + str(group_no) + ".faa"
        database_file = os.path.join(db_dir, db_name)
        with open(database_file, 'w') as fasta_out:
            for assembly in assembly_list:
                with gzip.open(os.path.join(assembly_dir,assembly+".faa.gz"), "rt") as infile:
                    shutil.copyfileobj(infile, fasta_out)
        db_size = "The estimated size of database is " + str(
            stat(database_file).st_size
        )
        cluster_report = os.path.join(db_dir, "cluster_report.txt")
        with open(cluster_report, "a", encoding="utf-8") as fin:
            fin.write("\n")
            fin.write(db_size + "\n")

def uniq_proteins(d_dir: str, db_name: str):
    """
    generate the protein search database with unique protein sequences
    :param d_dir:  path of directory containing databases
    :param db_name: name of the search database
    """
    unique_records = {}
    protein_file = os.path.join(d_dir, db_name)
    uniq_db = os.path.join(d_dir, "unique_" + db_name)
    for record in SeqIO.parse(protein_file, "fasta"):
        # if str(record.id) not in unique_records.keys():
        if str(record.id) not in unique_records:
            unique_records[str(record.id)] = record.seq
    with open(uniq_db, "w") as fout:
        for k, v in unique_records.items():
            fout.write(">" + k + "\n")
            fout.write(str(v) + "\n")

def remove_file(file_name):
    try:
        os.remove(file_name)
    except OSError as e:  ## if failed, report it back to the user ##
        print ("Error: %s - %s." % (e.filename, e.strerror))
