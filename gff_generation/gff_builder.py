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
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

def protein_report_processing(protein_report, protein_list: list,pride_id: str):
    """
    processing of peptide reports to yield a dataframe
    :param pride_id: ID of the associated metaproteomics study
    :param protein_list: list of expressed proteins
    :param protein_report: processed protein report
    """
    all_info=[]
    pep_info=pd.DataFrame(columns=['digest','contig_name','protein_start','protein_end','strand','Attributes','Validated PSMs' ,'Spectrum Counting'])
    for item in protein_list:
        item = str(item)
        all_info=[]
        temp_data=(protein_report.loc[protein_report['digest']==item]).reset_index(drop=True)
        unambiguous_peptides=[]
        ambiguous_peptides=[]
        vpsm=[]
        sc=[]
        Unique_peptide_to_protein_mapping="None"
        Ambiguous_peptide_to_protein_mapping="None"
        for i in range(len(temp_data)):
            if (int(temp_data["#Proteins"][i] == 1) and int(temp_data["Validated Protein Groups"][i]==1)):
                unambiguous_peptides.append(temp_data["Sequence"][i])
                vpsm.append(temp_data["Validated PSMs"][i])
                sc.append(temp_data["Spectrum Counting"][i])
                contig=temp_data["contig_name"][i]
                start=temp_data["protein_start"][i]
                end=temp_data["protein_end"][i]
                strand=temp_data["strand"][i]
            else:
                ambiguous_peptides.append(temp_data["Sequence"][i])
                vpsm.append(temp_data["Validated PSMs"][i])
                sc.append(temp_data["Spectrum Counting"][i])
                contig=temp_data["contig_name"][i]
                start=temp_data["protein_start"][i]
                end=temp_data["protein_end"][i]
                strand=temp_data["strand"][i]
        unambiguous_peptides = list(set(unambiguous_peptides))
        ambiguous_peptides = list(set(ambiguous_peptides))
        if len(unambiguous_peptides)==0 and len(ambiguous_peptides)>=1:
            all_info.append( (str(item) ,contig,start,end,strand, "ID="+str(item)+";type=Protein;Unique_peptide_to_protein_mapping=None;Ambiguous_peptide_to_protein_mapping=True;ambiguous_sequences="+",".join(ambiguous_peptides)+";pride_id="+pride_id,vpsm,sc) )
        elif len(unambiguous_peptides)>=1 and len(ambiguous_peptides)==0:
            all_info.append( (str(item) ,contig,start,end,strand, "ID="+str(item)+";type=Protein;Unique_peptide_to_protein_mapping=True;unambiguous_sequences="+",".join(unambiguous_peptides)+";Ambiguous_peptide_to_protein_mapping=None;pride_id="+pride_id,vpsm,sc) )
        else:
            all_info.append( (str(item) ,contig,start,end,strand, "ID="+str(item)+";type=Protein;Unique_peptide_to_protein_mapping=True;unambiguous_sequences="+",".join(unambiguous_peptides)+";Ambiguous_peptide_to_protein_mapping=True;ambiguous_sequences="+",".join(ambiguous_peptides)+";pride_id="+pride_id,vpsm,sc) )
        df = pd.DataFrame(all_info, columns=["digest", 'contig_name','protein_start','protein_end','strand','Attributes','Validated PSMs' ,'Spectrum Counting'])
        pep_info = pd.concat([pep_info,df], ignore_index=True)
    return pep_info

def gff_generation(attributes_file: str, assembly_name:str, out_folder: str):
    gff_data=attributes_file[['contig_name','protein_start','protein_end','strand','Attributes']]
    gff_data['strand'] = gff_data['strand'].map({-1: '-', 1: '+'})
    gff_data = gff_data.rename(columns={'contig_name':'seqid','protein_start':'start','protein_end':'end','strand':'strand','Attributes':'attributes'})
    gff_data = gff_data.assign(source='PeptideShaker')
    gff_data = gff_data.assign(type='CDS')
    gff_data = gff_data.assign(score='.')
    gff_data = gff_data.assign(phase='.')
    cols = gff_data.columns.to_list()
    col=cols[:1] + cols[5:7] +cols[1:3] +cols[7:8] +cols[3:4] +cols[-1:] +cols[4:5]
    gff_data=gff_data[col]
    topline='##gff-version 3'
    df_flatten = list(zip(*map(gff_data.get, gff_data)))
    out_file=os.path.join(out_folder, assembly_name+".gff")
    with open(out_file,'w', buffering=1) as out_handle:
        print('##gff-version 3', file=out_handle)
        for row in df_flatten:
            print('\t'.join([str(val) for val in row]), file=out_handle)
