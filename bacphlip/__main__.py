from .bacphlip import *

if __name__ == '__main__':
    args = parse_cmd_line_args()
    run_pipeline(args.input_file, force_overwrite=args.force_overwrite, local_hmmsearch=args.local_hmmsearch)


