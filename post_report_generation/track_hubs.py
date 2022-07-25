import pandas as pd


def get_sequence_header(input_fasta):
    header = []

    with open(input_fasta, 'r') as reader:
        original_proteins = reader.readlines()
        for line in original_proteins:
            if line[0] == '>':
                header.append(line[1:-1])

    return header


def remove_human_crap_protein_groups(proteinGroup):
    PGRemovalHumanAndCrap = ""
    tmp_pg = proteinGroup.replace(' ', '').split(',')
    for protein in tmp_pg:
        if not protein.startswith('MGYP9'):
            PGRemovalHumanAndCrap += protein
            PGRemovalHumanAndCrap += ","

    if PGRemovalHumanAndCrap.endswith(','):
        PGRemovalHumanAndCrap = PGRemovalHumanAndCrap[0:-1]

    return PGRemovalHumanAndCrap


def get_first_value(value_list):
    first_value = value_list[0]

    return first_value


def remove_doubtful_protein_groups(proteinGroups):
    PGLists = proteinGroups.replace(' ', '').split(';')
    PGAfterfilter = ""
    for pg in PGLists:
        if pg.endswith('(Confident)'):
            tmp_pg = pg[0:-11].replace(' ', '').split(',')
            PGRemovalHumanAndCrap = ""
            for protein in tmp_pg:
                if not protein.startswith('MGYP9'):
                    PGRemovalHumanAndCrap += protein
                    PGRemovalHumanAndCrap += ","

            if PGRemovalHumanAndCrap.endswith(','):
                PGRemovalHumanAndCrap = PGRemovalHumanAndCrap[0:-1]

            PGAfterfilter += PGRemovalHumanAndCrap
            PGAfterfilter += ";"

    if PGAfterfilter.endswith(';'):
        PGAfterfilter = PGAfterfilter[0:-1]

    return PGAfterfilter


def count_validated_protein_groups(processedPG):
    PGCount = 0
    if len(processedPG) > 0:
        PGCount = len(processedPG.split(';'))

    return PGCount


def remove_irrelevant_proteins(proteins, positions, proteinGroups):
    PG_proteins = []
    for pg in proteinGroups:
        for p in pg.replace(' ', '').split(','):
            PG_proteins.append(p)
    PG_proteins = list(set(PG_proteins))

    temp_proteinIDs = []
    temp_positions = []

    for i in range(len(proteins)):
        if proteins[i] in PG_proteins:
            temp_proteinIDs.append(proteins[i])
            temp_positions.append(positions[i])

    return temp_proteinIDs, temp_positions


def get_spectrum_counting(proteinGroups, SC_df):
    PGLists = proteinGroups.replace(' ', '').split(';')
    SC = ""
    for pg in PGLists:
        if pg in SC_df.keys():
            SC += str(round(SC_df[pg]['Spectrum Counting'],2))
            SC += ";"

    if SC.endswith(';'):
        SC = SC[0:-1]

    return SC


def create_dict_for_spectrum_counting(protein_report):
    protein = pd.read_csv(protein_report, sep='\t')[['Protein Group','Spectrum Counting']]
    protein['Processed Protein Group'] = protein['Protein Group'].apply(remove_human_crap_protein_groups)
    df = protein[['Processed Protein Group', 'Spectrum Counting']]
    df = df.drop_duplicates(subset='Processed Protein Group', keep="last")
    df = df.set_index('Processed Protein Group').to_dict(orient='index')

    return df


def apply_SC(peptides, SC_df):
    return pd.Series([get_spectrum_counting(ppg, SC_df) for ppg in peptides['Processed Protein Groups']])


def get_track_beds(peptide_report, protein_report, save_file_name, pxd_id):
    peptides = pd.read_csv(peptide_report, sep='\t')[['Protein(s)','Protein Group(s)','Sequence','Position','#Validated PSMs']]
    peptides = peptides[peptides['#Validated PSMs'] > 0]

    peptides['Processed Protein Groups'] = peptides['Protein Group(s)'].apply(remove_doubtful_protein_groups)
    peptides['Validated Protein Groups'] = peptides['Processed Protein Groups'].apply(count_validated_protein_groups)
    peptides = peptides[peptides['Validated Protein Groups'] > 0]
    peptides = peptides.reset_index()

    SC_df = create_dict_for_spectrum_counting(protein_report)
    # peptides['Spectrum Counting'] = peptides['Processed Protein Groups'].apply(get_spectrum_counting)
    peptides['Spectrum Counting'] = apply_SC(peptides, SC_df)

    proteins = peptides['Protein(s)'].str.replace(' ', '').str.split(';')
    proteinGroups = peptides['Processed Protein Groups'].str.replace(' ', '').str.split(';')
    positions = peptides['Position'].str.replace(' ', '').str.split(';')
    sequences = peptides['Sequence']
    validated_psms = peptides['#Validated PSMs']
    validated_PG = peptides['Validated Protein Groups']
    processed_PG = peptides['Processed Protein Groups']
    processed_SC = peptides['Spectrum Counting']


    processed_proteins = []
    processed_positions = []
    for i in range(len(proteins)):
        temp_proteinIDs, temp_positions = remove_irrelevant_proteins(proteins[i], positions[i], proteinGroups[i])
        processed_proteins.append(temp_proteinIDs)
        processed_positions.append(temp_positions)

    output_proteins = []
    output_sequences = []
    output_positions = []
    output_validated_psms = []
    output_validated_PG = []
    output_processed_PG = []
    output_processed_SC = []
    for i in range(len(processed_proteins)):
        protein_list = processed_proteins[i]
        seq = sequences[i]
        pos_list = processed_positions[i]
        valid_psm = validated_psms[i]
        valid_PG = validated_PG[i]
        pg = processed_PG[i]
        sc = processed_SC[i]

        for j in range(len(protein_list)):
            output_proteins.append(protein_list[j])
            output_sequences.append(seq)
            output_positions.append(pos_list[j])
            output_validated_psms.append(valid_psm)
            output_validated_PG.append(valid_PG)
            output_processed_PG.append(pg)
            output_processed_SC.append(sc)


    output = pd.DataFrame()
    output['Protein'] = output_proteins
    output['Sequence'] = output_sequences
    output['Position'] = output_positions
    output['Validated PSMs'] = output_validated_psms
    output['Validated Protein Groups'] = output_validated_PG
    output['Processed Protein Groups'] = output_processed_PG
    output['Spectrum Counting'] = output_processed_SC

    seq_count = []
    for i in range(len(output_proteins)):
        count = 0
        for j in range(i, len(output_proteins)):
            if output_proteins[i] == output_proteins[j] and output_sequences[i] == output_sequences[j]:
                count += 1
        seq_count.append(count)

    output['Protein Sequence Counts'] = seq_count
    output = output[output['Protein Sequence Counts'] == 1]

    groupby_sequence = output.groupby(['Sequence']).count()['Protein']
    seq_count = pd.DataFrame()
    seq_count['Sequence'] = groupby_sequence.index.to_list()
    seq_count['#Proteins'] = groupby_sequence.values
    output = output.merge(seq_count, how='left', on='Sequence')

    output['PXD ID'] = pxd_id
    output['PRIDE Link'] = 'ebi.ac.uk/pride/archive/projects/' + pxd_id

    save_columns = ['Protein','Sequence','Position','#Proteins','Validated Protein Groups','Validated PSMs','Spectrum Counting','Processed Protein Groups','PXD ID','PRIDE Link']
    output[save_columns].to_csv(save_file_name, index=False)
