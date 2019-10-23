#!/usr/bin/env python3

import os, sys
from ont_fast5_api import __version__
from argparse import ArgumentParser
from ont_fast5_api.multi_fast5 import MultiFast5File
from ont_fast5_api.conversion_tools.conversion_utils import get_fast5_file_list
   
def batch_reverter(input_path, output_folder, filename_base, batch_size, threads, recursive,
                   keys=set(('Raw', 'channel_id', 'context_tags', 'tracking_id'))):
    # make sure output dir doesn't exists
    if os.path.exists(output_folder):
        sys.stderr.write("Directory exists: %s\n"%output_folder)
        sys.exit(1)
    os.makedirs(output_folder)
    # get files to process - in revert order, since fail is typically before pass
    file_list = get_fast5_file_list(input_path, recursive)
    file_list = file_list[::-1]
    print("%s files to process..."%len(file_list))
    fi, ri = 0, -1
    for i, input_file in enumerate(file_list, 1):
        with MultiFast5File(input_file, 'r') as input_f5:
            for ri, read in enumerate(input_f5.get_read_ids(), ri+1):
                sys.stderr.write(" %s %s %s %s  \r"%(fi, ri, read, input_file))
                if not ri%batch_size:
                    output_f5 = MultiFast5File(os.path.join(output_folder, "%s_%s.fast5"%(filename_base, fi)), 'w')
                    fi += 1
                # copy group to new file
                read_name = "read_" + read
                group = input_f5.handle[read_name]
                output_f5.handle.copy(group, read_name)
                # and remove additional info
                reverted_group = output_f5.handle[read_name]#; print(reverted_group.keys())
                for k in reverted_group.keys():
                    if k not in keys:
                        del reverted_group[k]

def main():
    parser = ArgumentParser("Tool for reverting reads from a multi_read_fast5_file to raw information only")
    parser.add_argument('-i', '--input_path', required=True,
                        help='Folder containing single read fast5 files')
    parser.add_argument('-s', '--save_path', required=True,
                        help="Folder to output MultiRead subset to")
    parser.add_argument('-f', '--filename_base', default='batch', required=False,
                        help="Root of output filename, default='batch' -> 'batch_0.fast5'")
    parser.add_argument('-n', '--batch_size', type=int, default=4000, required=False,
                        help="Number of reads per multi-read file")
    parser.add_argument('-r', '--recursive', action='store_true',
                        help="Search recursively through folders for for single_read fast5 files")
    parser.add_argument('-t', '--threads', type=int, default=1, required=False,
                        help="Number of threads to use")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    args = parser.parse_args()
    
    batch_reverter(args.input_path, args.save_path, args.filename_base, args.batch_size, args.threads, args.recursive)

    
if __name__ == '__main__':
    main()
