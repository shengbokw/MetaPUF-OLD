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

import gzip
import os
import logging
import shutil
import mg_toolkit
import subprocess
import sys
import sourmash
import screed
import time
from argparse import ArgumentParser
from collections import defaultdict
import pandas as pd

from metagenomics_db import fetch_data as fd
from metagenomics_db import get_clusters as gc


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main():  # noqa: C901
    """
    Generates protein databases containing non-redundant proteins from
    metagenomics and/or metatranscriptomics assemblies
    """
    parser = ArgumentParser(
        description="Generate protein sequence database for Metagenomics and /or Metatranscriptomics datasets"
    )
    parser.add_argument(
        "-s","--study",
        type=str,
        required=True,
        help="Secondary study accession of the assembled study starting with ERP/SRP/DRP ",
    )
    parser.add_argument(
        "-v","--ver",
        type=str,
        required=True,
        help="pipeline version of MGnify analysis for the study",
    )
    parser.add_argument(
        "-i","--input_dir",
        type=dir_path,
        required=True,
        help="full path of the study folder containing raw folder",
    )
    parser.add_argument(
        "-m","--metadata",
        type=str, required=True,
        help="full path of the sample-assembly mapping file (.tsv) with filename."
    )

    starttime = time.time()
    args = parser.parse_args()
    sample_assembly_map = defaultdict(list)
    os.chdir(args.input_dir)
    fd.check_study_accession(args.study)
    samples = args.metadata
    with open(samples, "r") as f:
        headers = (f.readline()).strip().split("\t")
        sample_index = headers.index("sample_alias")
        analysis_index = headers.index("analysis_accession")
        for line in f:
            str1 = line.strip().split("\t")
            sample_assembly_map[str1[sample_index]].append(str1[analysis_index])
    #getting the sequnece data and predicted cds
    cmd_get_data = "  ".join(["mg-toolkit -d bulk_download -a",  args.study, "-p ", args.ver,  "-g sequence_data"])
    print(cmd_get_data)
    subprocess.call(cmd_get_data, shell=True)
    assembly_folder = os.path.join(args.input_dir, "assemblies")
    if not os.path.isdir(assembly_folder):
        subprocess.Popen(" ".join(["mkdir ", assembly_folder]), shell=True)
    os.makedirs(assembly_folder, exist_ok=True)
    sequence_dir=args.input_dir+"/"+args.study+"/5.0/sequence_data"
    for k, v in sample_assembly_map.items():
        sample_file = os.path.join(assembly_folder, k + ".fasta.gz")
        print(sample_file)
        with gzip.open(sample_file, "wt") as wf:
            for assembly in v:
                #renaming headers of the predicted cds in assemblies
                fd.rename_contigs(assembly, assembly_folder,sequence_dir)
                with gzip.open(os.path.join(sequence_dir,assembly+"_FASTA.fasta.gz"), "rt") as infile:
                    shutil.copyfileobj(infile, wf)
    matrix_file = os.path.join(assembly_folder, args.study+".tsv")

    meta_genomes = [assembly_folder+"/"+k+".fasta.gz" for k in sample_assembly_map.keys()]
    minhashes = []
    for contigs in meta_genomes:
        mh = sourmash.MinHash(ksize=31, n=0, scaled=1000)
        for record in screed.open(contigs):
            query_seq=record.sequence
            mh.add_sequence(query_seq, True)
        minhashes.append(mh)
    with open(matrix_file, 'w') as f_out:
        for i, e in enumerate(minhashes):
            basename=os.path.basename(meta_genomes[i])
            Name=".".join((basename).split(".")[:-2])
            _ = f_out.write(Name+',')
            for j, e2 in enumerate(minhashes):
                x = e.jaccard(minhashes[j])
                _ = f_out.write(str(round(x, 3))+',')
            _= f_out.write('\n')
    data = pd.read_csv(matrix_file, sep=",",header=None, index_col=[0])
    data=data.dropna(axis='columns', how='all')
    data.index.names=['index']
    col_names = data.index.to_list()
    data.columns=col_names
    # returns a dictionary containing assembly name and count of proteins in the asssembly
    
    proteins_info = fd.count_proteins(sample_assembly_map,sequence_dir)
    database_folder = os.path.join(assembly_folder, "databases")
    if not os.path.isdir(database_folder):
        subprocess.Popen(" ".join(["mkdir ", database_folder]), shell=True)
    os.makedirs(database_folder, exist_ok=True)
    # returns a dictionary with the group number and assemblies in the group
    samples_in_cluster = gc.generate_clusters(data, proteins_info,args.study, database_folder, sample_assembly_map)
    logging.info(f"Samples in the cluster: f{samples_in_cluster}")
    fd.build_db(args.study, database_folder,assembly_folder, samples_in_cluster)
    for file in os.listdir(database_folder):
        if file.endswith(".faa") and not file.startswith("unique"):
            protein_file=os.path.join(database_folder, file)
            fd.uniq_proteins(database_folder, file)
            fd.remove_file(protein_file)
    logging.info("Completed")
    logging.info("Runtime is {} seconds".format(time.time() - starttime))
if __name__ == "__main__":
    log_file = "db_generate.log"
    logging.basicConfig(
        level=logging.DEBUG, filemode="w", format="%(message)s", datefmt="%H:%M:%S"
    )
    main()

