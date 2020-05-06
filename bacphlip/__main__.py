from .bacphlip import *

if __name__ == '__main__':
    import argparse 
    ###Command line arguments    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",\
            required=True, help="Should be a valid path to a single genome (nucleotide) FASTA file containing only 1 record/contig.")
    parser.add_argument("-f", "--force_overwrite", action="store_true",\
            help="Whether to overwrite all existing files that will be created if they exist. Default is False")
    parser.add_argument("--local_hmmsearch", default=False,\
            help="By default, BACPHLIP assumes a system install of \"hmmsearch\". Use this flag to provide a custom path "
                    "to a local install of hmmsearch if necessary.")
    args = parser.parse_args()    
    run_pipeline(args.input_file, force_overwrite=args.force_overwrite, local_hmmsearch=args.local_hmmsearch)


