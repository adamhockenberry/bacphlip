#!/usr/bin/env python

import bacphlip
import argparse 

def main():
    args = bacphlip.parse_cmd_line_args()
    print(args)
    if args.multi_fasta:
        bacphlip.run_pipeline_multi(args.input_file, force_overwrite=args.force_overwrite, local_hmmsearch=args.local_hmmsearch)
    else:
        bacphlip.run_pipeline(args.input_file, force_overwrite=args.force_overwrite, local_hmmsearch=args.local_hmmsearch)

