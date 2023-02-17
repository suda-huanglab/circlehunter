from argparse import ArgumentParser
import textwrap
import sys
import os


import numpy as np
import pandas as pd


def run(
    mock_regions, out, size=1000, seed=0, multiple=6,
    length_loc=12, length_scale=3.5,
    length_minimum=5000, length_maximum=5_000_000,
    segments_loc=1, segments_scale=1
):
    # let all regions in a single axis
    mock_regions = pd.read_table(
        mock_regions, names=['chrom', 'start', 'end'], usecols=[0, 1, 2]
    )
    mock_regions['length'] = mock_regions['end'] - mock_regions['start']
    mock_regions['last'] = mock_regions['length'].cumsum()
    mock_regions['first'] = mock_regions['last'] - mock_regions['length']

    # make random regions
    rs = np.random.RandomState(seed)
    random_first = rs.randint(0, mock_regions.iloc[-1]['last'], multiple * size)
    rs = np.random.RandomState(seed)
    random_length = rs.normal(length_loc, length_scale, multiple * size)
    random_length = np.round(np.exp(random_length)).astype(int) + length_minimum
    random_length = random_length % length_maximum
    rs = np.random.RandomState(seed)
    random_strand = rs.choice(['+', '-'], multiple * size)
    rs = np.random.RandomState(seed)
    random_segments = rs.normal(segments_loc, segments_scale, multiple * size)
    random_segments = np.round(np.exp(random_segments)).astype(int) + 1
    no = np.repeat(
        np.arange(1, random_segments.shape[0] + 1), random_segments
    )[:multiple * size]
    random_region = pd.DataFrame({
        'first': random_first, 'length': random_length, 'strand': random_strand, 'no': no
    })

    # filter out regions fall in two clean region
    start_matrix = np.less.outer(
        random_region['first'].values, mock_regions['last'].values
    )
    random_region['start_index'] = np.argmax(start_matrix, axis=1)
    end_matrix = np.less.outer(
        random_region['first'].values + random_region['length'].values,
        mock_regions['last'].values
    )
    random_region['end_index'] = np.argmax(end_matrix, axis=1)
    random_region = random_region[
        random_region['start_index'] == random_region['end_index']
    ].copy()
    # projecting back to chromosome axis
    random_region = random_region.join(
        mock_regions[['chrom', 'first', 'start']], on='start_index', rsuffix='_ref'
    )
    random_region['start'] += random_region['first'] - random_region['first_ref']
    random_region['end'] = random_region['start'] + random_region['length']
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
    parser.add_argument('<mock_regions>', help='mock ecDNA fall in regions')
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
        '--length-maximum', dest='length_maximum', help='length maximum', default=5_000_000, type=int
    )
    parser.add_argument(
        '--segments-loc', dest='segments_loc', help='ln segments loc', default=1, type=float
    )
    parser.add_argument(
        '--segments-scale', dest='segments_scale', help='ln segments scale', default=1, type=float
    )
    args = vars(parser.parse_args())
    run(
        args['<mock_regions>'], args['<out>'],
        size=args['size'], seed=args['seed'], multiple=args['multiple'],
        length_loc=args['length_loc'], length_scale=args['length_scale'],
        length_minimum=args['length_minimum'], length_maximum=args['length_maximum'],
        segments_loc=args['segments_loc'], segments_scale=args['segments_scale']
    )


if __name__ == '__main__':
    main()
