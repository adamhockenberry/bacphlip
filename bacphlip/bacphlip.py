import pkg_resources
from Bio import SeqIO
import os
#import pandas as pd
#import joblib

SKLEARN_CLASSIFIER = pkg_resources.resource_filename('bacphlip', 'data/rf_best.joblib')
HMMER_DB = pkg_resources.resource_filename('bacphlip', 'data/prot_models.hmm')

def check_existing_file(infile):
    if not os.path.exists(infile):
        print("Input file {} does not appear to exist. Exiting.".format(infile))
        exit(1)

def check_nonexisting_file(outfile):
    if os.path.exists(outfile):
        print("Specified output file ({}) appears to already exist. Remove before running program again. Exiting.".format(args.output_file))
        exit(1)

def six_frame_translate(fasta_file_path, output_file_path, min_prot_length=40):
    """
    Should ensure that input is nucleotide data for six frame translation to make sense.
    """
    ###Basic error testing of input/output files
    check_existing_file(fasta_file_path)
    check_nonexisting_file(output_file_path)
    ###More error testing to ensure that the input file contains (the correct) data
    nt_records = list(SeqIO.parse(args.input_file, 'fasta'))
    if len(nt_records) == 1:
        nt_record = nt_records[0]
    elif len(nt_records) == 0:
        print('Input fasta file appears to be empty. Exiting.')
        exit(1)
    else:
        print('Input fasta file appears to contain more than one sequence record. Currently bacphlip only supports single contig inputs. Exiting.')
        exit(1)
    
    ###Run basic code
    genome_id = nt_record.id
    prots = []
    for i in [0, 1, 2]:###Each of the three reading frames
        ###Translate the regular strand
        tempy = nt_record[i:].translate()
        ###Split according to stop codons
        seq = str(tempy.seq).split('*')
        ###Append the sequences if they are longer than the defined minimum
        for j in seq:
            if len(j) >= min_prot_length:
                prots.append(j)
        
        ###Repeat for the reverse complement
        tempy = nt_record.reverse_complement()[i:].translate()
        seq = str(tempy.seq).split('*')
        for j in seq:
            if len(j) >= min_prot_length:
                prots.append(j)
    
    ###Write the results
    with open(output_file_path, 'w') as outfile:
        for i, seq in enumerate(prots):
            outfile.write('>{}_{}\n{}\n'.format(genome_id, i, seq))
    return    

#def process_args():
#    
#    return
#
#def check_args():
#    
#    return

#def predict_lifestyle(classifier, df):
#    ###Load classifier model
#    clf = joblib.load('../Data/rf_highMinAJH.joblib')
#    ###Load dataset
#    single_df = pd.read_csv(args.input_file, sep='\t', index_col=0)
#    ###Predict
#    class_probs = clf.predict_proba(single_df)
#    ###Write output
#    with open(args.output_file, 'w') as outfile:
#        outfile.write('{}\t{}\n'.format('Lytic', 'Temperate'))
#        outfile.write('{}\t{}\n'.format(class_probs[0][0], class_probs[0][1]))
#    exit(0)
#    return
#
#def process_single():
#
#    return
#
#def process_batch():
#
#    return

def main():
    import argparse 
    ###Command line arguments    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",\
            required=True, help="Should be a valid path to a single genome (nucleotide) fasta file containing only 1 record/contig.")
    parser.add_argument("-o", "--output_file",\
            required=True, help="Can be any valid path that does not currently exist.")
    args = parser.parse_args()
    six_frame_translate(args.input_file, args.output_file)
    exit(0)



if __name__ == '__main__':
    main()

#    ###
#    process_args()
#    
#    ###
#    check_args()
#    
#    for 
#    six_frame_translate()
#    hmmsearch_py()
#    process_hmmsearch()
#    predict_lifestyle()
#    pass
