import argparse
import sys
import pandas as pd
import subprocess


def get_args():
    """Collect the inputs for the executable scripts."""
    parser = argparse.ArgumentParser(
        description="""This script is going to map the samples, raw files and searching databases
            between PRIDE and MGnify, then exeucte shell scripts. """
    )

    parser.add_argument('-in', '--input', dest='input_file',
                        help="the sample information file, where stores the mapping information for samples, databases and raw files")
    parser.add_argument('-out', '--output', dest='output_folder',
                        help="the folder path for saving output files")
    parser.add_argument('-jar', '--javafile', dest='jar_file',
                        help="the mapping heads for the protein file")
    parser.add_argument('-par', '--parameters', dest='par_file',
                        help="the parameter file for searchgui script")
    args = parser.parse_args()

    if args.input_file is None:
        sys.exit('Error: no input file (sample information file)')
    if args.output_folder is None:
        sys.exit('Error: no output path provided!')
    if args.jar_file is None:
        sys.exit('Error: no jar file provided!')
    if args.par_file is None:
        sys.exit('Error: parameter file is mandatory!')

    return args


# input_file = 'config/sample_info.csv'
# par_file = 'config/rawurls.txt'
def searchgui_search(jar_file, input_file, output_folder, par_file):
    """
    run searchgui shell script within a loop based on the mapping sample information
    :input_file:    sample information file contains raw files, samples, mapped searching databases
    :jar_file:      the searchgui jar file
    :output_folder: output folder path for searchgui search results
    :par_file:      the searchgui parameters file
    """
    sample_info = pd.read_csv(input_file, sep=',')
    samples = sample_info['Sample'].drop_duplicates().to_list()

    for sample in samples:
        params = " -spectrum_files"
        params += " input/Raw/" + sample
        params += " -fasta_file"
        fasta = sample_info.loc[sample_info['Sample'] == sample, 'Database'].iloc[0]
        params += " assembilies/" + fasta + "_concatenated_target_decoy.fasta"
        params += " -output_folder"
        params += " " + output_folder
        params += " -id_params"
        params += " " + par_file
        params += " -xtandem 1 -msgf 1"
        params += " -output_default_name"
        params += " " + sample + "_searchgui"
        params += " -output_option 0 -output_data 1 -output_date 0"
        params += " -log logs/SearchGUI_search &> "
        params += "logs/" + sample + "_SearchGUI_search.log"
        params += " && touch "
        params += "searchgui/" + sample + "_searchgui.zip"

        commandline = "java -cp " + jar_file + " eu.isas.searchgui.cmd.SearchCLI "
        commandline += params

        subprocess.run(commandline)


def main():
    """
    Main function
    """
    args = get_args()
    searchgui_search(args.jar_file, args.input_file, args.output_folder, args.par_file)


if __name__ == "__main__":
    main()
java -cp /Users/shengbo_wang/Documents/MetaPUF/MetaPUF/workflow/bin/SearchGUI-4.0.41/SearchGUI-4.0.41.jar eu.isas.searchgui.cmd.SearchCLI  -spectrum_files input/Raw/S03 -fasta_file assembilies/ERZ4291133.fasta -output_folder searchgui -id_params searchgui/PXD020692_searchgui.par -xtandem 1 -msgf 1 -output_default_name S03_searchgui -output_option 0 -output_data 1 -output_date 0 -log logs/SearchGUI_search &> log/S03_SearchGUI_search.log && touch searchgui/S03_searchgui.zip
