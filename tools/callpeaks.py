# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: callpeaks.py
# time: 2021/02/16
from argparse import ArgumentParser


def read_bdg(bdg):
    with open(bdg) as f:
        for line in f:
            chrom, start, end, value = line.strip().split('\t')
            start, end, value = int(start), int(end), float(value)
            yield chrom, start, end, value


def bdg_filter(regions, cutoff):
    for chrom, start, end, value in regions:
        if value >= cutoff:
            yield chrom, start, end


def peak_merge(regions, gap):
    last_chrom, last_start, last_end = next(regions)
    for chrom, start, end in regions:
        if last_chrom == chrom and start - last_end <= gap:
            last_end = end
        else:
            yield last_chrom, last_start, last_end
            last_chrom, last_start, last_end = chrom, start, end
    else:
        yield last_chrom, last_start, last_end


def peak_filter(regions, length):
    for chrom, start, end in regions:
        if end - start >= length:
            yield chrom, start, end


def call_peaks(bdg, cutoff, gap, length, peak):
    regions = read_bdg(bdg)
    regions = bdg_filter(regions, cutoff)
    regions = peak_merge(regions, gap)
    regions = peak_filter(regions, length)
    with open(peak, 'w') as f:
        for chrom, start, end in regions:
            f.write(f'{chrom}\t{start}\t{end}\n')


def main():
    parser = ArgumentParser(
        description='Call peaks from bdg file'
    )
    parser.add_argument(
        '-c', dest='cutoff', help='cutoff of bases become a peak', default=2, type=float
    )
    parser.add_argument(
        '-g', dest='gap', help='maxmum gap between significant bases in a peak', default=200, type=int
    )
    parser.add_argument(
        '-l', dest='length', help='minimum length of peak', default=1000, type=int
    )
    parser.add_argument('<bdg>', help='input bdg file')
    parser.add_argument('<peak>', help='output peaks in bed format')
    args = vars(parser.parse_args())
    call_peaks(args['<bdg>'], args['cutoff'], args['gap'],
               args['length'], args['<peak>'])


if __name__ == '__main__':
    main()
