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


import json
import logging
import os

import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage


def generate_clusters(  # noqa: C901
    df1, prot_dict, study_name, p_dir, sample_assembly_info
):
    """
    generate cluster of assemblies to generte the db
    :param df1: dataframe containing sourmash comparison matrix
    :param prot_dict: dictinary with assembly name and protein counts
    :param study_name: name of the study (ERP/SRP/DRP)
    :param p_dir: directory path for cluster_report
    :param sample_assembly_info: dictionary containing sample assembly information
    """
    #df1=pd.read_csv(df1, sep=',', index_col=[0])
    col_name = df1.columns.to_list()
    #df1 = df1.set_axis(col_name, axis="index")
    group_no = 1
    max_clusters = len(prot_dict)
    assemblies = col_name
    cluster_assemblies = {}
    cluster_elements = {}
    cluster_report = os.path.join(p_dir, "cluster_report.txt")
    for cluster in range(1, max_clusters):
        if len(df1) < 1:
            break
        else:
            Z = linkage(df1, "ward")
            label = fcluster(Z, cluster, criterion="maxclust")
            df_clst = pd.DataFrame()
            df_clst["index"] = df1.index
            df_clst["label"] = label
            for i in range(cluster):
                elements = df_clst[df_clst["label"] == i + 1]["index"].tolist()
                size = len(elements)
                # size_of_cluster=str(size)+"_"+str(i+1)
                logging.info(f"\n Cluster {(i+1)}: N = {size}  {elements}")
                protein_counts = [prot_dict[x] for x in elements]
                group_size = 185 * sum(i for i in protein_counts)
                if group_size > 1073741824:
                    continue  # taking to the outer loop
                if group_size <= 1073741824:
                    logging.info(
                        f"cluster {cluster} found {size} assemblies for database {group_no}"
                    )
                    cluster_elements[group_no] = elements
                    group_no += 1
                    remaining_elements = list(
                        set(sorted(assemblies)).difference(set(sorted(elements)))
                    )
                    for item in elements:
                        temp_df = df1.loc[df1.index != item]
                        df1 = temp_df.copy()
                    new_data = df1[remaining_elements]
                    df1 = new_data.copy()
                    logging.info(len(df1))
                    assemblies = remaining_elements
    for k, v in cluster_elements.items():
        temp_list = []
        for item in v:
            # if item in sample_assembly_info.keys():
            if item in sample_assembly_info:
                for val in sample_assembly_info[item]:
                    temp_list.append(val)
        cluster_assemblies[k] = temp_list
    study_info = ("cluster_report for " + str(study_name)).center(40)
    database_info = "No of databases:" + str(group_no - 1)
    with open(cluster_report, "w") as fin:
        fin.write(study_info + "\n")
        fin.write(database_info + "\n")
        fin.write(json.dumps(cluster_assemblies))
    return cluster_assemblies
