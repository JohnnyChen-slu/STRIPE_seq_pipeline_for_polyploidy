#!/usr/bin/env python

import argparse
import gzip
import time
from Bio import SeqIO

def fastq_parser(input_file):
    with gzip.open(input_file, 'rt') as f:
        for record in SeqIO.parse(f, 'fastq'):
            yield record

def generate_filename(input_file, tag):
    base_name = input_file.split('.fastq.gz')[0]
    if "_R1_" in base_name:
        base_name = base_name.replace("_R1_", "_")
    elif "_R2_" in base_name:
        base_name = base_name.replace("_R2_", "_")
    return base_name + tag + '.fastq.gz'


def parse_args():
    parser = argparse.ArgumentParser(description='Process some fastq files.')
    parser.add_argument('input_file_r1', type=str, help='input file R1')
    parser.add_argument('input_file_r2', type=str, help='input file R2')
    return parser.parse_args()

def main():
    args = parse_args()

    tags = ['_good_reads_R1', '_good_reads_R2', '_bad_umi_R1', '_bad_umi_R2', '_bad_reads_R1', '_bad_reads_R2']
    files = [generate_filename(args.input_file_r1, tag) for tag in tags]

    good_reads = 0
    bad_reads = 0
    bad_umi = 0

    print(f"Processing the {args.input_file_r1.split('.fastq.gz')[0]} file...")

    start_time = time.time()

    with gzip.open(files[0], 'wt') as good_r1, \
         gzip.open(files[1], 'wt') as good_r2, \
         gzip.open(files[2], 'wt') as bad_umi_r1, \
         gzip.open(files[3], 'wt') as bad_umi_r2, \
         gzip.open(files[4], 'wt') as bad_r1, \
         gzip.open(files[5], 'wt') as bad_r2:

        r1_parser = fastq_parser(args.input_file_r1)
        r2_parser = fastq_parser(args.input_file_r2)

        for r1, r2 in zip(r1_parser, r2_parser):
            if str(r1.seq)[8:15] == 'TATAGGG':
                good_reads += 1
                SeqIO.write(r1, good_r1, 'fastq')
                SeqIO.write(r2, good_r2, 'fastq')
            elif 'TATAGGG' in str(r1.seq)[:20]:
                bad_umi += 1
                SeqIO.write(r1, bad_umi_r1, 'fastq')
                SeqIO.write(r2, bad_umi_r2, 'fastq')
            else:
                bad_reads += 1
                SeqIO.write(r1, bad_r1, 'fastq')
                SeqIO.write(r2, bad_r2, 'fastq')

    end_time = time.time()
    elapsed_time = end_time - start_time

    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    
    print('Finished!')
    print(f"Total good reads, Total bad reads, Total bad UMIs, elapsed time")
    print(f"{good_reads}, {bad_reads}, {bad_umi}, {int(hours)}:{int(minutes)}:{seconds:.2f}")

if __name__ == '__main__':
    main()

