from argparse import ArgumentParser
from multiprocessing import Pool
import textwrap
import sys
import os


from pybedtools import BedTool
from pyfaidx import Fasta
import numpy as np
import pandas as pd


def max_consecutive_n(ecDNA, fasta):
    """
    Calculate the maximum number of consecutive N

    Args:
        ecDNA (pd.Series): Series contain a region
        fasta (str): Path to an indexed fasta file
    """    
    def consecutive_n(fasta, chrom, start, end):
        """
        Find all consecutive N

        Args:
            fasta (str): Path to an indexed fasta file
            chrom (str): chromosome
            start (int): start
            end (end): end

        Yields:
            int: length of consecutive N
        """
        seq = Fasta(fasta, as_raw=True)[chrom][start:end]
        count = 0
        for base in seq:
            if base in ('n', 'N'):
                count += 1
            elif count > 0:
                yield count
                count = 0
        yield count

    return max(consecutive_n(fasta, ecDNA['chrom'], ecDNA['start'], ecDNA['end']))


def get_seq(genome, chrom, start, end, strand):
    genome = Fasta(genome)
    if strand == '+':
        return genome[chrom][start:end].seq
    elif strand == '-':
        return (-genome[chrom][start:end]).seq
    else:
        raise ValueError(f'{strand} is no validate argument for strand')


def run(
    chrom_sizes, blacklist, fasta, out, size=1000, seed=0, multiple=6,
    length_loc=12, length_scale=3.5, length_minimum=5000,
    segments_loc=1, segments_scale=1, max_n=10_000,
    processes=1
):
    # extract regions from the genome but exclude blacklist
    blacklist = BedTool(blacklist)
    clean_bed = blacklist.complement(g=chrom_sizes)
    clean_region = clean_bed.sort(g=chrom_sizes).to_dataframe()
    clean_region = clean_region[
        clean_region['chrom'].str.match(r'^chr[0-9X]{1,2}$')
    ].copy().reset_index(drop=True)
    # let all regions in a single axis
    clean_region['length'] = clean_region['end'] - clean_region['start']
    clean_region['last'] = clean_region['length'].cumsum()
    clean_region['first'] = clean_region['last'] - clean_region['length']
    # make random regions
    rs = np.random.RandomState(seed)
    random_first = rs.randint(
        0, clean_region.iloc[-1]['last'], multiple * size)
    random_length = np.round(np.exp(rs.normal(
        length_loc, length_scale, multiple * size))).astype(int) + length_minimum
    random_strand = rs.choice(['+', '-'], multiple * size)
    segments = np.round(np.exp(
        rs.normal(segments_loc, segments_scale, multiple * size))).astype(int) + 1
    no = np.repeat(
        np.arange(1, segments.shape[0] + 1), segments)[:multiple * size]
    random_region = pd.DataFrame({
        'first': random_first, 'length': random_length, 'strand': random_strand, 'no': no
    })
    # filter out regions fall in two clean region
    start_matrix = np.less.outer(
        random_region['first'].values, clean_region['last'].values
    )
    random_region['start_index'] = np.argmax(start_matrix, axis=1)
    end_matrix = np.less.outer(
        random_region['first'].values + random_region['length'].values,
        clean_region['last'].values
    )
    random_region['end_index'] = np.argmax(end_matrix, axis=1)
    random_region = random_region[random_region['start_index']
                                  == random_region['end_index']].copy()
    # projecting back to chromosome axis
    random_region = random_region.join(
        clean_region[['chrom', 'first', 'start']], on='start_index', rsuffix='_ref'
    )
    random_region['start'] += random_region['first'] - \
        random_region['first_ref']
    random_region['end'] = random_region['start'] + random_region['length']
    # count N length
    with Pool(processes=processes) as pool:
        n = pool.starmap(
            max_consecutive_n,
            ((ecDNA, fasta) for _, ecDNA in random_region.iterrows()),
            chunksize=100
        )
        random_region['consecutive_n'] = n
    # filter out region with too many N
    random_region = random_region[
        (random_region['consecutive_n'] < max_n)
        & (random_region['consecutive_n'] < random_region['length'])
    ].copy()
    mock_size = random_region['no'].value_counts().shape[0]
    if mock_size < size:
        print(
            f'Generated only {mock_size} ecDNA, you may need to increase "multiple"!',
            file=sys.stderr
        )
    # name each ecDNA
    random_region['no'] = random_region.groupby('no').ngroup() + 1
    random_region = random_region.join(
        random_region['no'].value_counts().rename('total'), on='no')
    random_region['part'] = random_region.groupby('no').cumcount() + 1
    random_region['name'] = (
        'ecDNA_' + random_region['no'].astype(str) + '_' +
        random_region['total'].astype(
            str) + '_' + random_region['part'].astype(str)
    )
    # save mock ecDNA
    real = random_region[random_region['no'] <= size].copy()
    real[['chrom', 'start', 'end', 'name', 'length', 'strand']].to_csv(
        f'{out}.bed', sep='\t', header=False, index=False
    )
    # generate bed for extract fasta
    if not os.path.isdir(f'{out}.fa') or not os.path.exists(f'{out}.fa'):
        os.mkdir(f'{out}.fa')
    for no, ecDNA in real.groupby('no'):
        cut = ecDNA.iloc[0]['length'] // 2
        first = ecDNA.iloc[0].copy()
        full = ecDNA.iloc[1:]
        last = ecDNA.iloc[0].copy()
        # shift the first segment to generate junction sequence
        if ecDNA.iloc[0]['strand'] == '+':
            first['start'] += cut + 500
            last['end'] -= cut
        elif ecDNA.iloc[0]['strand'] == '-':
            first['end'] -= cut - 500
            last['start'] += cut
        else:
            raise ValueError(
                f'{ecDNA.iloc[0]["strand"]} is no validate argument for strand')
        # save bed files
        bed = pd.concat([first.to_frame().T, full,
                         last.to_frame().T], ignore_index=True)
        bed[['chrom', 'start', 'end', 'name', 'length', 'strand']].to_csv(
            f'{out}.fa/ecDNA_{no}.bed', sep='\t', header=False, index=False
        )


def main():
    parser = ArgumentParser(description='mock ecDNA')
    parser.add_argument('<chrom_sizes>', help='input bam files')
    parser.add_argument('<blacklist>', help='output bed file')
    parser.add_argument('<fasta>', help='output bam file')
    parser.add_argument('<out>', help='output bam file')
    parser.add_argument(
        '-n', dest='size', help='how many ecDNA to mock', default=1000, type=int
    )
    parser.add_argument(
        '-m', dest='multiple', help='multiple', default=6, type=int
    )
    parser.add_argument(
        '-s', dest='seed', help='random seed', default=None, type=int
    )
    parser.add_argument(
        '-p', dest='processes', help='processes', default=1, type=int
    )
    parser.add_argument(
        '--length-loc', dest='length_loc', help='ln length loc', default=12, type=float
    )
    parser.add_argument(
        '--length-scale', dest='length_scale', help='ln length scale', default=3.5, type=float
    )
    parser.add_argument(
        '--length-minimum', dest='length_minimum', help='length minimum', default=5000, type=int
    )
    parser.add_argument(
        '--segments-loc', dest='segments_loc', help='ln segments loc', default=1, type=float
    )
    parser.add_argument(
        '--segments-scale', dest='segments_scale', help='ln segments scale', default=1, type=float
    )
    parser.add_argument(
        '--max-consecutive-n', dest='max_consecutive_n', help='max consecutive n', default=10_000, type=int
    )
    args = vars(parser.parse_args())
    run(
        args['<chrom_sizes>'], args['<blacklist>'], args['<fasta>'], args['<out>'],
        size=args['size'], seed=args['size'], multiple=args['multiple'],
        length_loc=args['length_loc'], length_scale=args['length_scale'], length_minimum=args['length_minimum'],
        segments_loc=args['segments_loc'], segments_scale=args['segments_scale'],
        max_n=args['max_consecutive_n'], processes=args['processes']
    )


if __name__ == '__main__':
    main()
