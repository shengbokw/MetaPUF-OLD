import os

##################################################
# SHELL
##################################################
# default executable for snakmake
shell.executable("bash")

##################################################
# PATHS
##################################################
# default configuration file
configfile:
    srcdir("config/config.proteomics.yaml")

# relevant paths
BINDIR      = srcdir("workflow/bin")
ENVDIR      = srcdir("workflow/envs")
OUTPUTDIR   = os.environ.get("OUTPUTDIR", config['outputdir'])
TMPDIR      = os.environ.get("TMPDIR", config['tmp_dir'])

##################################################
# WORKDIR
##################################################
workdir:
    OUTPUTDIR

#input files
I_PROTEINS = os.environ.get("PROTEINS", config["raws"]["Proteins"])
I_THERMORAW = os.environ.get("THERMORAW", config["raws"]["ThermoRaw"])
THERMOFOLD = os.environ.get("THERMOFOLD", config["raws"]["ThermoFold"])

# ##################################################
# # HARDWARE
# ##################################################
# BIGCORENO = os.environ.get("BIGCORENO", config["mem"]["big_mem_per_core_gb"])
# BIGMEMTOTAL = os.environ.get("BIGMEMTOTAL", config["mem"]["big_mem_total_gb"])
# MEMCORE = os.environ.get("MEMCORE", config["mem"]["normal_mem_per_core_gb"])
# BIGMEMCORE = os.environ.get("BIGMEMCORE", config["mem"]["big_mem_cores"])

# data
# PROTEINS_ORI = expand("proteins/{pname}_original.fasta", pname=I_PROTEINS)
# MAPPED_HEADER = expand("proteins/{pname}_mapping.csv", pname=I_PROTEINS)
PROTEINS_PROC  = expand("proteins/{pname}.fasta", pname=I_PROTEINS)
PROTEINS_DECOY = expand("proteins/{pname}_concatenated_target_decoy.fasta", pname=I_PROTEINS)

THERMORAW_NAMES = [os.path.splitext(os.path.basename(f))[0] for f in I_THERMORAW]
THERMOMGF = expand("{fname}/{bname}.mzML", fname=THERMOFOLD, bname=THERMORAW_NAMES)

SEARCHGUI_PAR  = expand("searchgui/{fname}_searchgui.par", fname=THERMOFOLD)
SEARCHGUI_ZIP  = expand("searchgui/{fname}_searchgui.zip", fname=THERMOFOLD)

# PEPTIDESHAKER_ZIP = expand("peptideshaker/{fname}_peptideshaker.cpsx", fname=THERMOFOLD)
# PEPTIDESHAKER_REPORT = expand("peptideshaker/{fname}_report.done", fname=THERMOFOLD)
PEPTIDESHAKER_MZID = expand("peptideshaker/{fname}_peptideshaker.mzid", fname=THERMOFOLD)
FINAL_MZID = expand("{fname}/{fname}_final.mzid", fname=THERMOFOLD)
PSM_TMP_RPT = expand("{fname}/peptideshaker_peptideshaker_1_Default_PSM_Report.txt", fname=THERMOFOLD)
PROTEIN_TMP_RPT = expand("{fname}/peptideshaker_peptideshaker_1_Default_Protein_Report.txt", fname=THERMOFOLD)
PEPTIDE_TMP_RPT = expand("{fname}/peptideshaker_peptideshaker_1_Default_Peptide_Report.txt", fname=THERMOFOLD)
PSM_RPT = expand("{fname}/{fname}_psm_report.txt", fname=THERMOFOLD)
PROTEIN_RPT = expand("{fname}/{fname}_protein_report.txt", fname=THERMOFOLD)
PEPTIDE_RPT = expand("{fname}/{fname}_peptide_report.txt", fname=THERMOFOLD)

# tools
THERMO_EXE = os.path.join(BINDIR, "ThermoRawFileParser/ThermoRawFileParser.exe")

SEARCHGUI_JAR = os.path.join(BINDIR, "SearchGUI-4.0.41/SearchGUI-4.0.41.jar")
SEARCHGUI_PAR_PARAMS = " ".join(["-%s %s" % (k, "'%s'" % v if isinstance(v, str) else str(v)) for k, v in config["searchgui"]["par"].items()])

PEPTIDESHAKER_JAR = os.path.join(BINDIR, "PeptideShaker-2.0.33/PeptideShaker-2.0.33.jar")

##################################################
# RULES
# Each subsequent output file needs to have its target path specified at the beginning.
##################################################
rule ALL:
    input:
        thermo=THERMOMGF,
        # mapping=[PROTEINS_PROC, MAPPED_HEADER],
        searchgui=[PROTEINS_DECOY, SEARCHGUI_PAR, SEARCHGUI_ZIP],
        report=[PROTEIN_TMP_RPT, PEPTIDE_TMP_RPT, PSM_TMP_RPT],
        peptideshaker=PEPTIDESHAKER_MZID
        # recover=[PROTEIN_RPT, PEPTIDE_RPT]


#########################
# ThermoRawFileParser
#########################
# https://github.com/compomics/ThermoRawFileParser
rule thermorawfileparser:
    input:
        exe=THERMO_EXE,
        raws=expand("{fname}/{bname}.raw",fname=THERMOFOLD,bname=THERMORAW_NAMES)
    output:
        THERMOMGF
    log:
        expand("logs/{fname}_thermorawfileparser.log",fname=THERMOFOLD)
    threads: 1
    message:
        "ThermoRawFileParser: {input} -> {output}"
    shell:
        "mono {input.exe} -d=$(dirname {input.raws[0]}) -o=$(dirname {output[0]}) -f=1 -m=0 &> {log}"


#########################
# Sequence header mapping change before SearchGUI
#########################
# rule replace_header:
#     input:
#         PROTEINS_ORI
#     output:
#         faa=PROTEINS_PROC,
#         map=MAPPED_HEADER
#     log:
#         expand("logs/{fname}_replace_header.log", fname=THERMOFOLD)
#     threads: 1
#     message:
#         "Replace header: {input} -> {output}"
#     shell:
#         "python replace_sequence_head.py -in {input} -out {output.faa} "
#         "-map {output.map} -b &> {log}"


#########################
# SearchGUI
#########################
rule searchgui_decoy:
    input:
        faa=PROTEINS_PROC,
        jar=SEARCHGUI_JAR
    output:
        PROTEINS_DECOY
    log:
        expand("logs/{fname}_SearchGUI_decoy.log",fname=THERMOFOLD)
    params:
        tmpdir = TMPDIR,
        logdir = "logs/SearchGUI_decoy"
    threads: 1
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    message:
        "SearchGUI decoy: {input} -> {output}"
    shell:
        "java -cp {input.jar} eu.isas.searchgui.cmd.FastaCLI -in {input.faa} "
        "-decoy -temp_folder {params.tmpdir} -log {params.logdir} &> {log}"


rule searchgui_config:
    input:
        jar=SEARCHGUI_JAR
    output:
        SEARCHGUI_PAR
    log:
        expand("logs/{fname}_SearchGUI_params.log",fname=THERMOFOLD)
    params:
        params = SEARCHGUI_PAR_PARAMS,
        tmpdir = TMPDIR,
        logdir = "logs/SearchGUI_params"
    threads: 1
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    message:
        "SearchGUI parameters: {input} -> {output}"
    shell:
        "java -cp {input.jar} eu.isas.searchgui.cmd.IdentificationParametersCLI -out {output} "
        "{params.params} -temp_folder {params.tmpdir} -log {params.logdir} &> {log}"


rule searchgui_search:
    input:
        par=SEARCHGUI_PAR,
        faa=PROTEINS_DECOY,
        mgf=THERMOMGF,
        jar=SEARCHGUI_JAR
    output:
        SEARCHGUI_ZIP
    log:
        expand("logs/{fname}_SearchGUI_search.log",fname=THERMOFOLD)
    params:
        name=expand("{fname}_searchgui", fname=THERMOFOLD),
        tmpdir = TMPDIR,
        logdir = "logs/SearchGUI_search"
    threads: 10
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    message:
        "SearchGUI search: {input.par}, {input.mgf} -> {output}"
    shell:
        """
        java -cp {input.jar} eu.isas.searchgui.cmd.SearchCLI \
            -spectrum_files $(dirname {input.mgf[0]}) \
            -fasta_file {input.faa} \
            -output_folder $(dirname {output}) \
            -id_params {input.par} \
            -xtandem 1 \
            -msgf 1 \
            -comet 0 \
            -andromeda 0 \
            -threads {threads} \
            -output_default_name {params.name} \
            -output_option 0 \
            -output_data 1 \
            -output_date 0 \
            -log {params.logdir} \
            &> {log} && touch {output}
        """


#########################
# PeptideShaker
#########################
# http://compomics.github.io/projects/peptide-shaker
rule peptideshaker_load:
    input:
        searchgui=SEARCHGUI_ZIP,
        jar=PEPTIDESHAKER_JAR
    output:
        protein=PROTEIN_TMP_RPT,
        peptide=PEPTIDE_TMP_RPT,
	psm=PSM_TMP_RPT,
        mzid=PEPTIDESHAKER_MZID
    log:
        expand("logs/{fname}_PeptideShaker_load.log",fname=THERMOFOLD)
    threads: 10
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    message:
        "PeptideShaker load SearchGUI results: {input.searchgui} -> {output.mzid}, {output.protein}, {output.peptide}, {output.psm}"
    shell:
        "java -cp {input.jar} eu.isas.peptideshaker.cmd.PeptideShakerCLI "
        "-reference 'peptideshaker_peptideshaker_1' "
        "-identification_files {input.searchgui} "
        "-out_reports $(dirname {output.protein}) -reports 3,6,9 "
        "-output_file {output.mzid} -contact_first_name 'Shengbo' -contact_last_name 'Wang' "
        "-contact_email 'shengbo_wang@ebi.ac.uk' -contact_address 'EBI' -organization_name 'EBI' "
        "-organization_email 'test@ebi.ac.uk' -organization_address 'Cambridge' "
        "-threads {threads} &> {log}"


#########################
# Recover sequence header sfter PeptideShaker
#########################
#rule recover_mzid:
#    input:
#        mzid=PEPTIDESHAKER_MZID,
#        map=MAPPED_HEADER
#    output:
#        FINAL_MZID
#    log:
#        expand("logs/{fname}_recover_header.log", fname=THERMOFOLD)
#    threads: 1
#    message:
#        "Replace header: {input} -> {output}"
#    shell:
#        "python replace_sequence_head.py -in {input.mzid} -out {output} "
#        "-map {input.map} -a &> {log}"


# rule recover_protein:
#     input:
#         map=MAPPED_HEADER,
#         protein=PROTEIN_TMP_RPT
#     output:
#         PROTEIN_RPT
#     log:
#         expand("logs/{fname}_recover_protein_report.log", fname=THERMOFOLD)
#     threads: 1
#     message:
#         "Replace header: {input.protein} -> {output}"
#     shell:
#         "python replace_sequence_head.py -in {input.protein} "
#         "-out {output} -map {input.map} -a &> {log}"


# rule recover_peptide:
#     input:
#         map=MAPPED_HEADER,
#         peptide=PEPTIDE_TMP_RPT
#     output:
#         PEPTIDE_RPT
#     log:
#         expand("logs/{fname}_recover_peptide_report.log", fname=THERMOFOLD)
#     threads: 1
#     message:
#         "Replace header: {input.peptide} -> {output}"
#     shell:
#         "python replace_sequence_head.py -in {input.peptide} "
#         "-out {output} -map {input.map} -a &> {log}"

# when using -zip to export files, the obtained zip file will not contain .cpsx
# file, what I did is first use -out to get the cpsx file, then use zip function
# to include the corresponding mgf files and fasta file.
# rule peptideshaker_report:
#     input:
#         archive=PEPTIDESHAKER_ZIP,
#         jar=PEPTIDESHAKER_JAR
#     output:
#         PEPTIDESHAKER_REPORT
#     log:
#         expand("logs/{fname}_PeptideShaker_report.log",fname=THERMOFOLD)
#     threads: 10
#     conda:
#         os.path.join(ENVDIR, "IMP_proteomics.yaml")
#     message:
#         "PeptideShaker report: {input.archive}"
#     shell:
#         "java -cp {input.jar} eu.isas.peptideshaker.cmd.ReportCLI "
#         "-in {input.archive} -out_reports $(dirname {output}) -reports 6,9 &> {log} && "
#         "touch {output}"
