# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: largeinsert.py
# time: 2021/02/16
from argparse import ArgumentParser
import pysam
import re


def read_chromosome_sizes(filename):
    chromosome_sizes = dict()
    with open(filename) as f:
        for line in f:
            chromosome, sizes = line.strip().split('\t')[:2]
            chromosome_sizes[chromosome] = int(sizes)
    return chromosome_sizes


def run(in_bam, out_bed, chromosome_sizes_file, include=2, exclude=1024, mapq=10, ratio=0.05, length=1500):
    chromosome_sizes = read_chromosome_sizes(chromosome_sizes_file)
    with pysam.AlignmentFile(in_bam, 'rb') as bam, open(out_bed, 'w') as bed:
        for chrom in chromosome_sizes.keys():
            for read in bam.fetch(chrom):
                if not read.flag & include:
                    continue
                if read.flag & exclude:
                    continue
                if read.mapq < mapq:
                    continue
                if read.has_tag('MD'):
                    mismatch = len(re.findall(r'[ATCG]', read.get_tag('MD')))
                    match = sum(
                        map(int, (re.findall(r'(\d+)M', read.cigarstring)))
                    )
                    if mismatch > match * ratio:
                        continue
                if read.reference_name != read.next_reference_name or abs(read.template_length) > length:
                    if read.is_reverse:
                        end = read.reference_end
                        start = max(0, end - length)
                        strand = '-'
                    else:
                        start = read.reference_start
                        end = min(
                            chromosome_sizes[read.reference_name], start + length
                        )
                        strand = '+'
                    bed.write(
                        f'{read.reference_name}\t{start}\t{end}\t{read.query_name}\t.\t{strand}\n'
                    )


def main():
    parser = ArgumentParser(
        description='Exract large insert site\'s reads span regions from Bam file.')
    parser.add_argument(
        '-l', dest='length', help='extend length of insert size', default=1500, type=int
    )
    parser.add_argument(
        '-f', dest='include', help='only include reads with all  of the FLAGs in INT present',
        default=1
    )
    parser.add_argument(
        '-F', dest='exclude', help='only include reads with none of the FLAGS in INT present',
        default=1028
    )
    parser.add_argument(
        '-q', dest='mapq', help='only include reads with mapping quality >= INT', default=10, type=int
    )
    parser.add_argument(
        '-r', dest='ratio', help='only include reads with mismatch ratio <= FLOAT', default=0.05, type=float
    )
    parser.add_argument('<chromosome_sizes>', help='chromosome sizes file')
    parser.add_argument('<in_bam>', help='input bam file')
    parser.add_argument('<out_bed>', help='output bed file')
    args = vars(parser.parse_args())
    run(
        args['<in_bam>'], args['<out_bed>'], args['<chromosome_sizes>'],
        args['include'], args['exclude'], args['mapq'], args['ratio'], args['length']
    )


if __name__ == '__main__':
    main()
