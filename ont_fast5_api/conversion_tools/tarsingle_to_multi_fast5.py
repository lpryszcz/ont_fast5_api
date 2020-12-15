#!/usr/bin/env python3
from argparse import ArgumentParser
from multiprocessing import Pool
from collections import deque
import logging
import os, sys, re

from ont_fast5_api import __version__
from ont_fast5_api.conversion_tools.conversion_utils import get_fast5_file_list, batcher, get_progress_bar
from ont_fast5_api.fast5_file import Fast5File
from ont_fast5_api.multi_fast5 import MultiFast5File
from single_to_multi_fast5 import add_read_to_multi_fast5

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
exc_info = False

import tarfile, urllib, urllib.request

# match all supported extensions
tarpat = re.compile('\.tar$|\.tgz$|\.tar\.gz$|\.tar\.gz.\d+$')

def get_progress_remote(f, filesize):
    # this results sometimes in AttributeError: '_io.BufferedReader' object has no attribute 'length'
    return 100.*(filesize-f.fileobj.fileobj.length)/filesize

def get_progress_local(f, filesize):
    return 100.*f.fileobj.fileobj.tell()/filesize    

def batch_convert_tarsingle_to_multi(input_path, save_path, filename_base, batch_size,
                                     tmp_dir, revert):
    """Convert archives with single Fast5 files into multi_fast5 files,
    optionally reverting to raw data (stripping all analyses from Fast5). 

    Note, for large tar archives, it'll occupy lots of memory with time. WHY??!!
    """
    # https://sra-pub-src-1.s3.amazonaws.com/SRR7415631/barcode05_bsubtilis_pass.tar.gz.1
    # ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR288/ERR2887851/Zymo-PromethION-LOG-BB-SN_signal.tar.gz
    # ~/cluster/dna_mods/ecoli/_archives/tar/$acc/barcode05_bsubtilis_pass.tar.gz.1
    m = tarpat.findall(os.path.basename(input_path))
    if not m:
        print("Unsupported file extension!")
        return
    # strip extension from file name
    output_folder = os.path.join(save_path, os.path.basename(input_path)[:-len(m[0])])
    # create outdir
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    else:
        print("Directory exists: %s"%output_folder)
        return
    print("Saving %s to %s ..."%(input_path, output_folder))
    # get mode
    mode = "r:gz" if input_path.endswith("gz") else "r"
    # open tar file, either remote or local, getting filesize
    if input_path.startswith(("ftp", "http", "www")):
        stream = urllib.request.urlopen(input_path)
        f = tarfile.open(fileobj=stream, mode=mode)
        filesize = int(stream.info()['Content-Length'])
        get_progress = get_progress_remote
    else:
        f = tarfile.open(input_path, mode)
        filesize = os.path.getsize(input_path)
        get_progress = get_progress_local
        
    fi = 0
    output_table = open(os.path.join(output_folder, "filename_mapping.txt"), 'w')
    output_table.write("single_read_file\tmulti_read_file\n")    
    for tarinfo in f:
        if os.path.splitext(tarinfo.name)[-1] == ".fast5":
            # report progress
            if not fi%100:
                sys.stderr.write(" %s reads \r"%fi)
                #sys.stderr.write(" %s reads [%5.1f%s] \r"%(fi, get_progress(f, filesize), '%'))
            # open new multi_fast5 file every batch_size of reads
            if not fi % batch_size:
                f.members = []
                multi_read_file = os.path.join(output_folder, "{}_{}.fast5".format(filename_base, int(fi/batch_size)))
                multi_f5 = MultiFast5File(multi_read_file, 'w')
            # extract to tmpdir, ideally to ramdrive so no disk I/O is involved
            single_read_file = tarinfo.name
            handle = os.path.join(tmp_dir, tarinfo.name)
            f.extract(tarinfo, path=tmp_dir)
            # get tar handle - reading from tar.gz is slooow, likely due to many seeks
            #handle = f.extractfile(tarinfo) #
            # open single & multi fast5
            try:
                single_f5 = Fast5File(handle, 'r')
            except Exception as e:
                print("[ERROR] at %s"%tarinfo.name)#, e)
                continue
            # store info
            add_read_to_multi_fast5(multi_f5, single_f5, revert)
            # store filename_mapping
            output_table.write("{}\t{}\n".format(single_read_file, os.path.basename(multi_read_file)))
            # clean up
            if os.path.isfile(handle):
                os.remove(handle)
            fi += 1
    output_table.close()
    print(" %s reads stored in %s Fast5 files.      "%(fi, int(fi/batch_size)+1))
            
def main():
    parser = ArgumentParser("")
    parser.add_argument('-i', '--input_path', nargs="+", 
                        help='TAR file containing single read fast5 files (can be gzipped and remote file ie ftp or http URL)')
    parser.add_argument('-s', '--save_path', required=True,
                        help="Folder to output multi read files to (input file basename will be added as well)")
    parser.add_argument('-f', '--filename_base', default='batch', required=False,
                        help="Root of output filename, default='batch' -> 'batch_0.fast5'")
    parser.add_argument('-n', '--batch_size', type=int, default=4000, required=False,
                        help="Number of reads per multi-read file")
    parser.add_argument('-t', '--tmp', default="/tmp", help="Folder from temporary files")
    parser.add_argument('--revert', action='store_true',
                        help="Revert Fast5 to original content. This will drop basecalling and other analyses.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    args = parser.parse_args()

    print("Processing %s archive(s)..."%len(args.input_path))
    for input_path in args.input_path:
          batch_convert_tarsingle_to_multi(input_path, args.save_path, args.filename_base,
                                           args.batch_size, args.tmp, args.revert)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
