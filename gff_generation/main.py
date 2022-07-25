# -*- coding: utf-8 -*-

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

import os
import shutil
import logging
import time
import subprocess
from argparse import ArgumentParser
import pandas as pd

from gff_generation import gff_builder as gb

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main():  # noqa: C901
    """
    Aggregate information from metaproteomics and metagenomics about the expressed proteins 
    and generate GFF format files for visualisation
    """
    parser = ArgumentParser(
        description="Aggregate information from metaproteomics and metagenomics about the expressed proteins and generate GFF format files for visualisation"
    )
    parser.add_argument(
        "-s","--sample_info",
        type=str,
        required=True,
        help="Absolute path of the sample metadata file",
    )
    parser.add_argument(
        "-r","--reports_dir",
        type=dir_path,
        required=True,
        help="full path of the directory containing processed metaproteomics reports",
    )
    parser.add_argument(
        "-m","--metag_dir",
        type=dir_path,
        required=True,
        help="full path of the directory containing contig metadata file",
    )
    parser.add_argument(
        "-p","--pride_id",
        type=str,
        required=True,
        help="PRIDE id for the reanalysed study",
    )

    starttime = time.time()
    args = parser.parse_args()
    sample_info=pd.read_csv(os.path.join(args.reports_dir,args.sample_info), sep=',')
    samples=list(set(sample_info['Sample'].to_list()))
    results_folder=os.path.join(args.reports_dir,"results")
    if not os.path.isdir(results_folder):
        subprocess.Popen(" ".join(["mkdir ", results_folder]), shell=True)

    os.makedirs(results_folder, exist_ok=True)
    sample_file_list = [args.reports_dir+"/"+f for f in os.listdir(args.reports_dir) if f.endswith('_peptide_report.csv')]
    csv_list = []
    for file in sorted(sample_file_list):
        csv_list.append(pd.read_csv(file))
    csv_merged = pd.concat(csv_list, ignore_index=True)
    csv_merged.to_csv(os.path.join(results_folder,"allsamples.csv"), index=False)
    proteins=list(set(csv_merged['Protein'].to_list()))
    unique_proteins=pd.DataFrame(proteins, columns=['digest'])
    assemblies=list(set(sample_info['Assembly'].to_list()))
    for assembly in assemblies:
        temp_assembly=pd.read_csv(os.path.join(args.metag_dir,assembly+"_contig_info.txt"), sep='\t')
        temp_assembly.columns=['assembly','contig_name','digest','partial_info', 'protein_start','protein_end','strand','caller']
        assembly_subset_of_proteins=pd.merge(unique_proteins, temp_assembly, on='digest', how='inner')
        assembly_expressed_proteins=assembly_subset_of_proteins.merge(csv_merged, left_on='digest', right_on='Protein',how='inner')
        assembly_expressed_proteins.to_csv(os.path.join(results_folder,assembly+"_expressed_proteins.csv"))
        expressed_proteins=list(set(assembly_expressed_proteins['digest']))
        attributes_file=gb.protein_report_processing(assembly_expressed_proteins,expressed_proteins,args.pride_id)
        gb.gff_generation(attributes_file, assembly, results_folder)

    logging.info("Completed")
    logging.info("Runtime is {} seconds".format(time.time() - starttime))


if __name__ == "__main__":
    log_file = "gff_generate.log"
    logging.basicConfig(
        level=logging.DEBUG, filemode="w", format="%(message)s", datefmt="%H:%M:%S"
    )
    main()
