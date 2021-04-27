# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: superfilter.py
# time: 2021/03/05
from argparse import ArgumentParser
import pandas as pd
import numpy as np

NAMES = [
    "chrom", "start", "end", "name", "score", "strand"
]


def run(typical, output):
    data = pd.read_table(typical, usecols=range(5), names=NAMES)
    data = data.sort_values('score', ascending=False)
    values = data['score'].values.copy()[::-1]
    values[values < 0] = 0
    minimum = np.sum(values <= np.quantile(values, 0.05))
    maximum = np.sum(values <= np.quantile(values, 0.95))
    slope = (values[maximum] - values[minimum]) / (maximum - minimum)
    x = np.arange(values.shape[0])
    b = values - (x * slope)
    y = np.add.outer(x * slope, b)
    count = np.sum(y.T >= values, axis=1)
    index = np.argmin(count)
    result = data[data['score'] >= values[index]]
    result.to_csv(output, sep='\t', index=False, header=False)


def main():
    parser = ArgumentParser(
        description='Filter bed with extreme high score.'
    )
    parser.add_argument('<input>', help='input bed file')
    parser.add_argument('<output>', help='output bed file')
    args = vars(parser.parse_args())
    run(args['<input>'], args['<output>'])


if __name__ == '__main__':
    main()
