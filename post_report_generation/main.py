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
import subprocess
from argparse import ArgumentParser
import pandas as pd

from post_report_generation import track_hubs as th

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main():  # noqa: C901
    """
    Aggregate information from metaproteomics and metagenomics about the expressed proteins
    and generate processed peptide reports for gff file generation
    """
    parser = ArgumentParser(
        description="Metadata for the post processing of peptide reports")
    parser.add_argument(
        "-s","--sample_info",
        type=str,
        required=True,
        help="Absolute path of the sample metadata file",)
    parser.add_argument(
        "-p","--pride_id",
        type=str,
        required=True,
        help="PRIDE id for the reanalysed study",)

    args = parser.parse_args()
    sample_info=pd.read_csv(args.sample_info, sep=',')
    samples=list(set(sample_info['Sample'].to_list()))

    for sample in samples:
        th.get_track_beds('results/reports/peptides/'+sample+'_peptide_report.txt',
                        'results/reports/proteins/'+sample+'_protein_report.txt',
                        'results/reports/processed/processed_'+sample+'_peptide_report.csv',
                        args.pride_id)


if __name__ == "__main__":
    main()
