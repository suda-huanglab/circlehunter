# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: bam2bed.py
# time: 2021/02/16
from argparse import ArgumentParser
import pysam
import re


def read_regions(filename):
    with open(filename) as f:
        for line in f:
            chrom, start, end = line.strip('\n').split('\t')[:3]
            start, end = int(start), int(end)
            yield chrom, start, end


def fetch_reads(bam, regions=None):
    if regions is None:
        yield from bam
    else:
        regions = read_regions(regions)
        for region in regions:
            for read in bam.fetch(*region):
                yield read


def run(in_bam, out_bed, regions=None, include=2, exclude=1024, mapq=10, ratio=0.05, insert_size=0):
    with pysam.AlignmentFile(in_bam, 'rb') as bam, open(out_bed, 'w') as bed:
        for read in fetch_reads(bam, regions):
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
            if read.reference_name == read.next_reference_name and abs(read.template_length) < insert_size:
                continue
            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end
            bed.write(f'{chrom}\t{start}\t{end}\n')


def main():
    parser = ArgumentParser(
        description='Extract read align positions as a BED format file.'
    )
    parser.add_argument(
        '-L', dest='regions', help='only include reads overlapping this BED FILE', default=None
    )
    parser.add_argument(
        '-f', dest='include', help='only include reads with all  of the FLAGs in INT present',
        default=1, type=int
    )
    parser.add_argument(
        '-F', dest='exclude', help='only include reads with none of the FLAGS in INT present',
        default=1028, type=int
    )
    parser.add_argument(
        '-q', dest='mapq', help='only include reads with mapping quality >= INT', default=10, type=int
    )
    parser.add_argument(
        '-r', dest='ratio', help='only include reads with mismatch ratio <= FLOAT', default=0.05, type=float
    )
    parser.add_argument(
        '-i', dest='insert_size', help='only include read insert size >= INT', default=0, type=int
    )
    parser.add_argument('<in_bam>', help='input bam file')
    parser.add_argument('<out_bed>', help='output bed file')
    args = vars(parser.parse_args())
    run(
        args['<in_bam>'], args['<out_bed>'], args['regions'], args['include'], args['exclude'], args['mapq'], args['ratio'], args['insert_size']
    )


if __name__ == '__main__':
    main()
