# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: circlehunter.py
# time: 2020/9/23
from collections import defaultdict
from itertools import groupby, count
from argparse import ArgumentParser
from pysam import AlignedSegment
from pysam import AlignmentFile
import re


M = 0
S = 4
H = 5

SPLIT_OPERATIONS = {S, H}


LEFT_CIGAR_REGEX = re.compile(r'\d+[HS]\d+M')

CENTER_CIGAR_REGEX = re.compile(r'\d+M')

RIGHT_CIGAR_REGEX = re.compile(r'\d+M\d+[HS]')

STRAND_DICT = {True: '-', False: '+'}


class InvalidReads(Exception):
    pass


class SplitReads:
    def __init__(self):
        self._left, self._center, self._right = None, None, None

    @property
    def is_full(self):
        return all(map(lambda x: x is not None, (self._left, self._center, self._right)))

    def add_read(self, read: AlignedSegment):
        """
        add the read into the SplitReads instance
        :param read: read
        :return: True if the SplitReads is full
        :raise InvalidReads: if the read does not meet the conditions
        """
        if self._left is None:
            self._validate_left(read)
            self._left = read
        elif self._center is None:
            self._validate_center(read)
            self._center = read
        elif self._right is None:
            self._validate_right(read)
            self._right = read
            self._left_shift()
        return self.is_full

    def _validate_left(self, read):  # noqa
        # left cigar must match nSmM
        if LEFT_CIGAR_REGEX.fullmatch(read.cigarstring) is None:
            raise InvalidReads(f'Invalid left cigar')
        # left has another one map position
        if not read.has_tag('SA') or read.get_tag('SA').count(';') != 1:
            raise InvalidReads(f'Reads has not another or more than two map positions')

    def _validate_center(self, read):
        # center map uniquely
        if read.has_tag('XA'):
            raise InvalidReads(f'Reads not mapped uniquely')
        if read.is_reverse == self._left.is_reverse:
            raise InvalidReads(f'Invalid center strand')
        # center full map(flag 0x2 ensures that mate is mapped and in different strand)
        # todo: support indel in center?
        if CENTER_CIGAR_REGEX.fullmatch(read.cigarstring) is None:
            raise InvalidReads(f'Invalid center cigar')
        # center is map righter
        if self._left.reference_start >= read.reference_start:
            raise InvalidReads(
                f'Invalid center reference start'
            )

    def _validate_right(self, read):
        # right cigar match nMmS
        if RIGHT_CIGAR_REGEX.fullmatch(read.cigarstring) is None:
            raise InvalidReads(f'Invalid right cigar')
        # right is in same strand
        if self._left.is_reverse != read.is_reverse:
            raise InvalidReads(f'Invalid right strand')
        # right map rightest
        if read.reference_start <= self._center.reference_start:
            raise InvalidReads(f'Invalid right map position')

    def _left_shift(self):
        new_tags = []
        align = self._left.query_alignment_length + self._right.query_alignment_length
        shift = align - self._left.infer_read_length()
        if shift > 0:
            new_tags.append(('LS', shift, 'i'))
            try:
                left_seq = self._left.query_alignment_sequence[:shift]
                right_seq = self._right.query_alignment_sequence[-shift:]
            except TypeError:
                print(self._right.query_name)
                return
            if left_seq == right_seq:
                new_tags.append(('MF', left_seq, 'Z'))
        self._right.set_tags(self._right.get_tags() + new_tags)

    def get_circle(self):
        if not self.is_full:
            raise RuntimeError(f'{self} is not full')
        shift = self._right.get_tag('LS') if self._right.has_tag('LS') else 0
        return (
            self._left.reference_name, self._left.reference_start, self._right.reference_end - shift
        )

    def add_circle_id(self, circle_id):
        for read in (self._left, self._center, self._right):
            new_tags = read.get_tags() + [('CI', circle_id, 'Z')]
            read.set_tags(new_tags)

    def get_reads(self):
        if not self.is_full:
            raise RuntimeError(f'{self} is not full')
        return self._left, self._center, self._right


def mapped_filter(reads):
    for read in reads:
        if read.flag & 2:
            yield read


def duplicate_filter(reads):
    for read in reads:
        if not read.is_duplicate:
            yield read


def fetch_split_reads(bam, contig):
    catches = {}

    # fetch reads from bam file
    reads = bam.fetch(contig)

    # filter out unmap and duplicate reads
    reads = mapped_filter(reads)
    reads = duplicate_filter(reads)

    for read in reads:
        # init or get a SplitReads instance
        if read.query_name not in catches:
            catches[read.query_name] = SplitReads()
        try:
            # try to add the read to the SplitReads instance
            full = catches[read.query_name].add_read(read)
            if full:  # split read is full, yield it
                yield catches.pop(read.query_name)
        except InvalidReads:  # invalid read occurred, clean up
            del catches[read.query_name]


def run(in_bam, out_circles, out_bam):
    bam = AlignmentFile(in_bam, 'rb')
    out = AlignmentFile(out_bam, 'wb', template=bam)

    circle_ids = count(1)
    circle_mapper = {}
    junction_counter = defaultdict(int)

    for contig in bam.references:
        for split_reads in fetch_split_reads(bam, contig):

            chrom, start, end = split_reads.get_circle()

            if (chrom, start, end) in circle_mapper:
                circle_id = circle_mapper[(chrom, start, end)]
            else:
                circle_id = f'C{next(circle_ids)}'
                circle_mapper[(chrom, start, end)] = circle_id
            junction_counter[circle_id] += 1
            split_reads.add_circle_id(circle_id)

            for read in split_reads.get_reads():
                out.write(read)

    bam.close()
    out.close()

    with open(out_circles, 'w') as f:
        for (chrom, start, end), circle_id in circle_mapper.items():
            f.write(f'{chrom}\t{start}\t{end}\t{circle_id}\t{junction_counter[circle_id]}\n')


def num_filter(reads):
    for _, reads_pair in groupby(reads, lambda read: read.query_name):
        reads_pair = tuple(reads_pair)
        if len(reads_pair) == 3:
            yield reads_pair


def chromosome_filter(reads):
    for pairs in reads:
        if len(set(read.reference_name for read in pairs)) == 1:
            yield pairs


def sort_pairs(reads):
    for pairs in reads:
        yield sorted(pairs, key=lambda read: read.reference_start)


def reads_pairs_filter(reads):
    for pairs in reads:
        if pairs[0].is_read1 != pairs[2].is_read1:
            continue
        if pairs[0].is_read1 == pairs[1].is_read1:
            continue
        yield pairs


def full_map_filter(reads):
    for pairs in reads:
        if len(pairs[1].cigartuples) == 1:
            yield pairs


def strand_filter(reads):
    for pairs in reads:
        if pairs[0].is_reverse == pairs[1].is_reverse:
            continue
        if pairs[0].is_reverse != pairs[2].is_reverse:
            continue
        yield pairs


def split_reads_filter(reads):
    for pairs in reads:
        if len(pairs[0].cigartuples) != 2 or pairs[0].cigartuples[1][0] != M:
            continue
        if len(pairs[2].cigartuples) != 2 or pairs[2].cigartuples[0][0] != M:
            continue
        yield pairs


def check_filter(reads):
    for pairs in reads:
        if pairs[0].cigartuples[0][1] < pairs[2].cigartuples[1][1]:
            yield pairs


def left_shift_reads(reads):
    for left, center, right in reads:
        shift = left.cigartuples[1][1] + right.cigartuples[0][1] - center.cigartuples[0][1]
        if shift > 0:
            left_seq = left.query_alignment_sequence[:shift]
            right_seq = right.query_alignment_sequence[-shift:]
            if left_seq == right_seq:
                right.set_tags(
                    right.get_tags() + [('LS', shift, 'i')]
                )
                yield left, center, right
        else:
            yield left, center, right


def main():
    parser = ArgumentParser(description='circlehunter')
    parser.add_argument('<in.bam>', help='input bam files')
    parser.add_argument('<out.bed>', help='output bed file')
    parser.add_argument('<out.bam>', help='output bam file')
    args = vars(parser.parse_args())
    run(args['<in.bam>'], args['<out.bed>'], args['<out.bam>'])


if __name__ == '__main__':
    main()
