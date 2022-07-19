# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: tssenrich.py
# time: 2020/11/18
from argparse import ArgumentParser
from multiprocessing import Pool
from functools import partial
import gzip

import numpy as np
import pysam


def get_coverage(
        bam_file, chrom, start, end, strand,
        bins=400, flag_require=2, flag_filter=1024, min_mapping_quality=20
):
    with pysam.AlignmentFile(bam_file) as bam:
        cols = bam.pileup(
            chrom, start, end, flag_require=flag_require,
            flag_filter=flag_filter, min_mapping_quality=min_mapping_quality
        )
        try:
            pos, depth = np.array([[col.pos, col.n] for col in cols]).T
        except ValueError:
            return np.zeros(bins)
    x = np.arange(start, end)
    y = np.zeros_like(x)
    indices = (pos >= start) & (pos < end)
    y[pos[indices] - start] = depth[indices]
    xi = np.linspace(start, end, bins)
    yi = np.interp(xi, x, y)
    return yi[::-1] if strand == '-' else yi


def run(
        bam_file, feature_file, extend=2000, bins=400, flag_require=2,
        flag_filter=1024, min_mapping_quality=20, processes=1
):
    target = partial(
        get_coverage, bam_file, bins=bins, flag_require=flag_require,
        flag_filter=flag_filter, min_mapping_quality=min_mapping_quality
    )
    pool = Pool(processes=processes)
    matrix = []
    # add features to queue
    if feature_file.endswith('.gz'):
        feature_reader = gzip.open(feature_file, 'rt')
    else:
        feature_reader = open(feature_file, 'rt')
    for line in feature_reader:
        if line.startswith('#'):
            continue
        chrom, start, end, _name, _score, strand = line.strip().split('\t')
        start, end = int(start) - extend, int(end) + extend
        pool.apply_async(target, args=(chrom, start, end, strand), callback=matrix.append)
    feature_reader.close()
    # wait for results
    pool.close()
    pool.join()
    matrix = np.vstack(matrix)
    # greenleaf style norm
    num_edge_bins = int(100 / (2 * extend / bins))
    bin_means = matrix.mean(axis=0)
    noises = np.sum(bin_means[:num_edge_bins]) + np.sum(bin_means[-num_edge_bins:])
    avg_noise = noises / (2 * num_edge_bins)
    matrix /= avg_noise
    tss_enrich = np.max(matrix.mean(axis=0))
    print(tss_enrich)


def main():
    parser = ArgumentParser(description='tssenrich')
    parser.add_argument('<in.bam>', help='input bam files')
    parser.add_argument('<feature>', help='feature bed file')
    parser.add_argument('<prefix>', help='output file prefix')
    parser.add_argument('-e', '--extend', help="extend length", default=2000)
    parser.add_argument('-b', '--bins', help="bins num", default=400)
    parser.add_argument(
        '-f', type=int, default=2, metavar='flag',
        help='only include reads with all  of the FLAGs in INT present'
    )
    parser.add_argument(
        '-F', type=int, default=1024, metavar='flag',
        help='only include reads with none of the FLAGS in INT present'
    )
    parser.add_argument(
        '-q', type=int, default=20, metavar='MAPQ',
        help='only include reads with mapping quality greater than or equal to'
    )
    parser.add_argument(
        '-p', type=int, default=1, metavar='processes', help='number of additional threads to use'
    )
    args = vars(parser.parse_args())
    run(
        args['<in.bam>'], args['<feature>'], extend=args['extend'], bins=args['bins'],
        flag_require=args['f'], flag_filter=args['F'], min_mapping_quality=args['q'],
        processes=args['p']
    )


if __name__ == '__main__':
    main()
