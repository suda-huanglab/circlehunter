# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: annotate.py
# time: 2020/10/21
from argparse import ArgumentParser
import pysam
import os


def fetch_names(db, chrom, start, end):
    try:
        return ','.join(record.split('\t')[3] for record in db.fetch(chrom, start, end))
    except ValueError:
        return ''


def main():
    parser = ArgumentParser(
        description='fetch database record and append names as the last column of input bed'
    )
    parser.add_argument('<db>', help='database')
    parser.add_argument('<input.bed>', help='input bed file')
    parser.add_argument('<output.bed>', help='output bed file')
    args = vars(parser.parse_args())

    db_path = os.path.abspath(args['<db>'])
    db = pysam.TabixFile(db_path)

    with open(args['<input.bed>']) as bed, open(args['<output.bed>'], 'w') as out:
        # header
        while True:
            line = bed.readline()
            if line.startswith('#'):
                out.write(line)
            else:
                break

        # append header
        out.write(f'## annotate database: {db_path}\n')

        # process records
        while line:
            line = line.strip()
            chrom, start, end, *_ = line.split('\t')
            names = fetch_names(db, chrom, int(start), int(end))
            out.write(f'{line}\t{names}\n')
            # read next line
            line = bed.readline()


if __name__ == '__main__':
    main()
