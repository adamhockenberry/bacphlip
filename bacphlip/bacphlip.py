import pkg_resources
from Bio import SeqIO, SearchIO
import os
import subprocess
import pandas as pd
from collections import OrderedDict
import joblib
import sys

"""
TODO:
    Check that input fastas are DNA and output is amino acids
    Check for suspiciously short/long genome lengths
    Allow for batch input (one fasta per line)
    Allow for existing amino acid file input (rather than 6frame)
    Allow for multi-contig input files
    Allow for custom hmmsearch
"""

SKLEARN_CLASSIFIER = pkg_resources.resource_filename('bacphlip', 'data/rf_best.joblib')
HMMER_DB = pkg_resources.resource_filename('bacphlip', 'data/prot_models.hmm')

TOTAL_PROTEIN_MODELS=206
MIN_PROT_LENGTH=40

def check_existing_file(infile):
    """
    Checks whether a file does not exist and provides a hopefully informative error message
    if that's not the case.

    Inputs:
        infile - path to a hopefully existing file
    """
    if not os.path.exists(infile):
        raise Exception("Input file ({}) does not appear to exist.".format(infile))

def check_nonexisting_file(outfile):
    """
    Checks whether a file to be written already  exists and provides an informative message if
    that check fails.
    
    Inputs:
        outfile - path to a hopefully non-existent file file
    """ 
    if os.path.exists(outfile):
        raise Exception("Specified output file ({}) appears to already exist. " 
                "Remove before running program again or specify force_overwrite flag ('-f').".format(outfile))

def check_genome_fasta_reqs(nt_records):
    """
    Checks whether a fasta file (meant to be a genome (i.e. DNA) file) can be read and ensures that it only
    contains one record.
    
    Another check to implement is that the fasta record is of DNA type.

    Inputs:
        nt_records - list from SeqIO.parse that ideally contains one element.
    """
    if len(nt_records) == 1:
        nt_record = nt_records[0]
    elif len(nt_records) == 0:
        raise Exception("Input fasta file appears to be empty. "
                "Ensure that the file contains what you think it does.")
    else:
        raise Exception("Input fasta file appears to contain more than one sequence record. "
                "Currently bacphlip only supports single contig inputs.")
    return nt_record

def six_frame_translate(fasta_file_path, output_file_path, force_overwrite=False):
    """
    Reads a genome fasta file and outputs an amino acid fasta file containing all possible
    six frame translational reading frames that are longer than a set length threshold.

    Inputs:
        fasta_file_path - system path to a valid genome fasta file
        output_file_path - system path for writing the amino acid file
        force_overwrite - if True will write "output_file_path" regardless of whether it exists
    """
    ###Basic error testing of input/output files
    check_existing_file(fasta_file_path)
    if not force_overwrite:
        check_nonexisting_file(output_file_path)
    ###Read in file    
    nt_records = list(SeqIO.parse(fasta_file_path, format='fasta'))
    nt_record = check_genome_fasta_reqs(nt_records)
    ###Run basic code
    genome_id = nt_record.id
    prots = []
    for i in [0, 1, 2]:###Each of the three reading frames
        modulo = len(nt_record[i:])%3
        assert modulo < 3
        ###Translate the regular strand
        if modulo > 0:
            tempy = nt_record[i:-modulo].translate()
        else:
            tempy = nt_record[i:].translate()
        ###Split according to stop codons
        seq = str(tempy.seq).split('*')
        ###Append the sequences if they are longer than the defined minimum
        for j in seq:
            if len(j) >= MIN_PROT_LENGTH:
                prots.append(j)
        ###Repeat for the reverse complement
        if modulo > 0:
            tempy = nt_record.reverse_complement()[i:-modulo].translate()
        else:
            tempy = nt_record.reverse_complement()[i:].translate()
        seq = str(tempy.seq).split('*')
        for j in seq:
            if len(j) >= MIN_PROT_LENGTH:
                prots.append(j)
    ###Write the results
    with open(output_file_path, 'w') as outfile:
        for i, seq in enumerate(prots):
            outfile.write('>{}_{}\n{}\n'.format(genome_id, i, seq))
    return    

def hmmsearch_py(aa_fasta_file, hmmsearch_out, force_overwrite=False, local_hmmsearch=False):
    """
    Runs hmmsearch on an amino acid fasta file against a pre-compiled database of protein domains
    of interest. 

    Inputs:
        aa_fasta_file - system path to a valid amino acid fasta file
        hmmsearch_out - system path for writing the hmmsearch output file
        force_overwrite - if True will write "hmmsearch_out" regardless of whether it exists
    """

    check_existing_file(aa_fasta_file)
    if not force_overwrite:
        check_nonexisting_file(hmmsearch_out)
    #
    hmmsearch_call = "hmmsearch"
    if local_hmmsearch:
        if not os.path.exists(local_hmmsearch):
            raise Exception("You specified a local hmmsearch install (path={}) but this file does not appear to exist.".format(local_hmmsearch))
        hmmsearch_call = local_hmmsearch
    #
    with open(hmmsearch_out, 'w') as outfile:
        try:
            subprocess.run(args=[hmmsearch_call, HMMER_DB, aa_fasta_file], check=True, stdout=outfile)
        except subprocess.CalledProcessError as e:
            raise Exception("System call to outside program \"hmmsearch\" failed. Please ensure that you have hmmsearch 3.x "\
            "or greater properly installed. BACPHLIP assumes by default that the install is available in the system path but "\
            "local installs can be supplied with the \"--local_hmmsearch\" flag followed by the path to a local hmmsearch.") from e
    return

def process_hmmsearch(hmmsearch_file, hmmsearch_df_out, force_overwrite=False):
    """
    Processes the hmmsearch default output file into a pandas dataframe -> tab-seperated file
    describing the presence/absence of each protein in the hmmsearch file.

    Inputs:
        hmmsearch_file - system path to a valid hmmsearch output file
        hmmsearch_df_out - system path for writing the processed tsv file
        force_overwrite - if True will write "hmmsearch_df_out" regardless of whether it exists
    """
    check_existing_file(hmmsearch_file)
    if not force_overwrite:
        check_nonexisting_file(hmmsearch_df_out)
    with open(hmmsearch_file, 'r') as infile:
        results = list(SearchIO.parse(infile, 'hmmer3-text'))
        simple_res = []
        for i in results:
            if len(i.hits) > 0:
                simple_res.append((i.id, 1))
            else:
                simple_res.append((i.id, 0))
    if len(simple_res) != TOTAL_PROTEIN_MODELS:
        raise Exception('Appears to be an error, too many or too few results given expected number of protein models tested.')
    simple_res = sorted(simple_res, key=lambda x: x[0])
    single_df = pd.DataFrame(OrderedDict(simple_res), index=[0])
    single_df.to_csv(hmmsearch_df_out, sep='\t')
    return

def predict_lifestyle(hmmsearch_df, predictions_out, force_overwrite=False):
    """
    Predicts the lifestyle of a phage given the location to a table describing the presence/absence of a set of
    protein domains and a pre-trained scikit-learn classifier.

    Inputs:
        hmmsearch_df - system path to a valid tsv file
        predictions_out - system path for writing the final predictions file
        force_overwrite - if True will write "hmmsearch_df_out" regardless of whether it exists
    """
    if not force_overwrite:
        check_nonexisting_file(predictions_out)
    clf = joblib.load(SKLEARN_CLASSIFIER)
    ###Load dataset
    single_df = pd.read_csv(hmmsearch_df, sep='\t', index_col=0)
    ###Predict
    class_probs = clf.predict_proba(single_df)
    ###Write output
    with open(predictions_out, 'w') as outfile:
        outfile.write('{}\t{}\n'.format('Virulent', 'Temperate'))
        outfile.write('{}\t{}\n'.format(class_probs[0][0], class_probs[0][1]))
    return



def run_pipeline(input_file_path, force_overwrite=False, local_hmmsearch=False):
    """
    Command-line implementation of the full BACPHLIP prediction pipeline. Currently implemented only for single genome
    inputs but lots of time could be saved during the classifier prediction step by implementing a batch option.
    """
    ### 
    six_frame_file = input_file_path + '.6frame'
    hmmsearch_file = input_file_path + '.hmmsearch'
    hmmsearch_df = input_file_path + '.hmmsearch.tsv'
    predictions_file = input_file_path + '.bacphlip'
    ###
    print('#################################################################################')
    print('Beginning BACPHLIP pipeline')
    six_frame_translate(input_file_path, six_frame_file, force_overwrite=force_overwrite)
    print('Finished six frame translation of genome (nucleotide) with output stored in {}'.format(six_frame_file))
    hmmsearch_py(six_frame_file, hmmsearch_file, force_overwrite=force_overwrite, local_hmmsearch=local_hmmsearch)
    print('Finished outside call to hmmsearch with output stored in {}'.format(hmmsearch_file))
    process_hmmsearch(hmmsearch_file, hmmsearch_df, force_overwrite=force_overwrite)
    print('Finished converting hmmsearch with output stored in {}'.format(hmmsearch_df))
    predict_lifestyle(hmmsearch_df, predictions_file, force_overwrite=force_overwrite)
    print('Finished with BACPHLIP predictions! Final output file stored in {}'.format(predictions_file))
    print('#################################################################################')

