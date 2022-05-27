# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: cutsites.py
# time: 2021/02/16
from argparse import ArgumentParser
import pysam

CUTSITE = {
    99: 'reference_start',
    147: 'reference_end',
    163: 'reference_start',
    83: 'reference_start'
}

SHIFT_FACTOR = {99: 4, 147: -5, 163: 4, 83: -5}


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


def cutsite(read):
    try:
        cutsite = getattr(read, CUTSITE[read.flag])
        cutsite += SHIFT_FACTOR[read.flag]
    except KeyError:
        return None
    return read.reference_name, cutsite


def run(in_bam, out_bed, regions=None, include=2, exclude=1024, mapq=10):
    with pysam.AlignmentFile(in_bam, 'rb') as bam, open(out_bed, 'w') as bed:
        for read in fetch_reads(bam, regions):
            if not read.flag & include:
                continue
            if read.flag & exclude:
                continue
            if read.mapq < mapq:
                continue
            record = cutsite(read)
            if record is not None:
                bed.write('{0}\t{1}\t{1}\n'.format(*record))


def main():
    parser = ArgumentParser(
        description='Extract Tn5 cutsite from mapped ATAC-Seq reads')
    parser.add_argument('-L',
                        dest='regions',
                        help='only include reads overlapping this BED FILE',
                        default=None)
    parser.add_argument(
        '-f',
        dest='include',
        help='only include reads with all  of the FLAGs in INT present',
        default=2,
        type=int)
    parser.add_argument(
        '-F',
        dest='exclude',
        help='only include reads with none of the FLAGS in INT present',
        default=1028,
        type=int)
    parser.add_argument('-q',
                        dest='mapq',
                        help='only include reads with mapping quality >= INT',
                        default=10,
                        type=int)
    parser.add_argument('<in_bam>', help='input bam file')
    parser.add_argument('<out_bed>', help='output bed file')
    args = vars(parser.parse_args())
    run(args['<in_bam>'], args['<out_bed>'], args['regions'], args['include'],
        args['exclude'], args['mapq'])


if __name__ == '__main__':
    main()
