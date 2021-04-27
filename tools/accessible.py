# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: accessible.py
# time: 2021/02/16
from argparse import ArgumentParser
import pysam

FRAGMENT = {
    99: 'reference_start', 147: 'reference_end', 163: 'reference_start', 83: 'reference_start'
}

SHIFT_FACTOR = {99: 4, 147: -5, 163: 4, 83: -5}


def read_chromosome_sizes(filename):
    chromosome_sizes = dict()
    with open(filename) as f:
        for line in f:
            chromosome, sizes = line.strip().split('\t')[:2]
            chromosome_sizes[chromosome] = int(sizes)
    return chromosome_sizes


def shift(read, chromosome_sizes, length=25):
    try:
        shift_pos = getattr(
            read, FRAGMENT[read.flag]) + SHIFT_FACTOR[read.flag]
    except KeyError:
        return None
    strand = '-' if read.is_reverse else '+'
    shift_start = max(shift_pos - length, 0)
    shift_end = min(shift_pos + length, chromosome_sizes[read.reference_name])
    return read.reference_name, shift_start, shift_end, '.', '.', strand


def run(in_bam, out_bed, chromosome_sizes_file, include=2, exclude=1024, mapq=10, length=25):
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
                record = shift(read, chromosome_sizes, length)
                if record is not None:
                    bed.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(*record))


def main():
    parser = ArgumentParser(
        description='Shift reads to cut sizes for ATAC-Seq')
    parser.add_argument(
        '-l', dest='length', help='extend length of cut sizes', default=25, type=int
    )
    parser.add_argument(
        '-f', dest='include', help='only include reads with all  of the FLAGs in INT present',
        default=2
    )
    parser.add_argument(
        '-F', dest='exclude', help='only include reads with none of the FLAGS in INT present',
        default=1028
    )
    parser.add_argument(
        '-q', dest='mapq', help='only include reads with mapping quality >= INT', default=10, type=int
    )
    parser.add_argument('<chromosome_sizes>', help='chromosome sizes file')
    parser.add_argument('<in_bam>', help='input bam file')
    parser.add_argument('<out_bed>', help='output bed file')
    args = vars(parser.parse_args())
    run(
        args['<in_bam>'], args['<out_bed>'], args['<chromosome_sizes>'],
        args['include'], args['exclude'], args['length']
    )


if __name__ == '__main__':
    main()
