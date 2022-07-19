# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: extractreads.py
# time: 2022/02/27
from argparse import ArgumentParser
import gzip


def run(fq_in, id_list, fq_out):
    fq_in = gzip.open(fq_in, 'rt')
    id_list = open(id_list, 'r')
    fq_out = gzip.open(fq_out, 'wt')

    next_id = id_list.readline().strip()
    while True:
        if not next_id:
            break
        reads = [fq_in.readline() for _ in range(4)]
        if not reads[0]:
            break
        reads_id = reads[0].split()[0]
        if reads_id == f'@{next_id}':
            fq_out.writelines(reads)
            next_id = id_list.readline().strip()

    fq_in.close()
    id_list.close()
    fq_out.close()


def main():
    parser = ArgumentParser(
        description='Extract reads from fastq file according to given id list')
    parser.add_argument('<fq_in>', help='input fastq.gz file')
    parser.add_argument('<id_list>', help='file contain reads id list')
    parser.add_argument('<fq_out>', help='output fastq.gz file')
    args = vars(parser.parse_args())
    run(args['<fq_in>'], args['<id_list>'], args['<fq_out>'])


if __name__ == '__main__':
    main()
