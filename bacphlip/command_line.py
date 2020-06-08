#!/usr/bin/env python

import bacphlip
import argparse 

def main():
    args = bacphlip.parse_cmd_line_args()
    bacphlip.run_pipeline(args.input_file, force_overwrite=args.force_overwrite, local_hmmsearch=args.local_hmmsearch)

