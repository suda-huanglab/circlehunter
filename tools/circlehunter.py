# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: circlehunter.py
# time: 2021/03/28
from collections import namedtuple, defaultdict, Counter, deque
from itertools import combinations
from functools import lru_cache
from argparse import ArgumentParser
import operator
import sys
import re

import numpy as np
import pandas as pd
import networkx as nx
import pysam


PEAK = 'PEAK'
LARGE = 'LARGE'

NEXT_CONNECTION_KEY = {
    PEAK: LARGE,
    LARGE: PEAK
}

SMALL_VALUE = np.exp(-10)
DIRECTION_INDEX = {'+': -1, '-': 0}
DIRECTION_FILTER = {'+': True, '-': False}
POS_ATTR = {'+': 'reference_end', '-': 'reference_start'}
OFFSET_FUNC = {'+': np.min, '-': np.max}
SPLIT_POS_FILTER = {'+': np.greater_equal, '-': np.less_equal}
HYPO_TRANS = {'+': np.add, '-': np.subtract}

Region = namedtuple('Region', ['chrom', 'start', 'end'])


def peaks2graph(peaks):
    """
    Generage a graph within the given largeinsert peaks.

    Args:
        peaks (string): path to the large peaks bed file.
    """
    peaks = pd.read_table(
        peaks, names=['chrom', 'start', 'end', 'peak', 'peak_weight', 'large_weight'], usecols=range(6)
    )

    graph = nx.MultiGraph()

    for _, *interval, peak, peak_weight, large_weight in peaks.itertuples():
        graph.add_node(
            Region(*interval), peak=peak, peak_weight=peak_weight, weight=large_weight
        )

    return graph


def count_reads_pairs(bam, regions, include=1, exclude=1036, mapq=10, ratio=0.05, length=1500):
    """
    count reads pairs span two different regions, which means these two region are connected

    Args:
        bam ([pysam.AlignmentFile]): pysam AlignmentFile object poited to the bam file
        regions (iterable): iterable object that contain regions that are large insert reads enrich
        include (int, optional): only include reads within these flags. Defaults to 1.
        exclude (int, optional): exclude reads within these flags. Defaults to 1036.
        mapq (int, optional): minimum mapq for accepted reads. Defaults to 10.
        ratio (float, optional): maximum ratio of mismatch. Defaults to 10.
        length (int, optional): minimum length for insert size that reads pair's insert size are large. Defaults to 1500.

    Returns:
        dict: link between two region and the supported read pairs count
    """
    # fetch reads from regions
    reads_pairs = defaultdict(lambda: defaultdict(list))
    for region in regions:
        for read in bam.fetch(*region):
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
            if read.reference_name != read.next_reference_name or abs(read.template_length) > length:
                pair = 1 if read.is_read1 else 2
                strand = '-' if read.is_reverse else '+'
                insert_size = read.template_length if read.template_length else sys.maxsize
                reads_pairs[read.query_name][pair].append(
                    (region, strand, insert_size)
                )

    # filter out reads is not paired and keep reads pair that is align far away
    reads_pairs = {
        name: {
            i: max(reads, key=lambda r: r[2])[:2] for i, reads in pair.items()
        }
        for name, pair in reads_pairs.items() if len(pair) > 1
    }

    # count links between regions by reads pair
    pairs_count = Counter(
        (tuple(pair.values()) for pair in reads_pairs.values())
    )

    return pairs_count.most_common()[::-1]


def add_large_connections(graph, bam, limit, include=1, exclude=1036, mapq=10, ratio=0.05, min_insert_size=1500):
    """
    add connection between two large insert enriched regions by paired reads

    Args:
        graph (networkx.MultiGraph): graph with large insert enriched regions as nodes
        bam ([pysam.AlignmentFile]): pysam AlignmentFile object poited to the bam file
        limit (number): the minimum supported paired reads to consider as a connection
        include (int, optional): only include reads within these flags. Defaults to 1.
        exclude (int, optional): exclude reads within these flags. Defaults to 1036.
        mapq (int, optional): minimum mapq for accepted reads. Defaults to 10.
        ratio (float, optional): maximum ratio of mismatch. Defaults to 10.
        min_insert_size (int, optional): minimum length for insert size that reads pair's insert size are large. Defaults to 1500.

    Returns:
        networkx.MultiGraph: graph that has been add new edges
    """
    # count reads pairs cross two enrich region
    connections = count_reads_pairs(
        bam, graph.nodes, include, exclude, mapq, ratio, min_insert_size
    )
    # connection nodes
    for ((node1, strand1), (node2, strand2)), count in connections:
        if count >= limit and node1 != node2:
            graph.add_edge(
                node1, node2, LARGE, weight=count, strand={
                    node1: strand1, node2: strand2
                }
            )

    return graph


def has_direction(graph, node, direction):
    """
    check the graph if the node has a large connection direction

    Args:
        graph (networkx.MultiGraph): graph with large insert enriched regions as nodes
        node (tuple): a large insert enrich region
        direction (str): direction

    Returns:
        bool: True if node has the direction connection
    """
    for _, edge in graph.adj[node].items():
        try:
            if edge[LARGE]['strand'][node] == direction:
                return True
        except KeyError:
            continue
    return False


def add_peak_connections(graph):
    """
    connect nodes within the same peak

    Args:
        graph (networkx.MultiGraph): graph with large insert enriched regions as nodes

    Returns:
        networkx.MultiGraph: graph has connection within the same peak
    """
    for (node1, meta1), (node2, meta2) in combinations(graph.nodes.items(), 2):
        if meta1['peak'] != meta2['peak']:
            continue
        # peak should match the true situation
        # the left side has a minus direction and right side should has a plus direction
        node1, node2 = sorted([node1, node2])
        if not has_direction(graph, node1, '-'):
            continue
        if not has_direction(graph, node2, '+'):
            continue
        graph.add_edge(
            node1, node2, PEAK, peak=meta1['peak'], weight=meta1['peak_weight']
        )
    return graph


def connect_types(graph, node):
    """
    find all connection types of the given node

    Args:
        graph (networkx.MultiGraph): graph with large insert enriched regions as nodes
        node (tuple): a large insert enrich region

    Yields:
        str: connection types
    """
    for _, meta in graph.adj[node].items():
        for key in meta.keys():
            yield key


def sorted_circle(nodes):
    """
    sort the nodes in a circle to make it unique

    Args:
        nodes (list): list of nodes that form a nodes

    Returns:
        tuple: tuple of nodes that form the input nodes but has been sorted
    """
    # the minimum node is the first node
    first = min(nodes)
    first_index = nodes.index(first)
    last_index = first_index - 1
    # reconstruct the cycle to make is fit the order but wont break the cycle
    next_index = first_index + 1 if first_index + 1 != len(nodes) else 0
    if nodes[last_index] > nodes[next_index]:
        nodes = nodes[first_index:] + nodes[:first_index]
    else:
        nodes = nodes[:next_index][::-1] + nodes[next_index:][::-1]
    return tuple(nodes)


def find_peak_child(graph, grandpa, father):
    """
    find next node which is connected by a peak connection

    Args:
        graph (networkx.MultiGraph): graph with large insert enriched regions as nodes
        grandpa (tuple): grandpa of the finding child
        father (tuple): father of the finding child

    Yields:
        tuple: children
    """
    validater = operator.lt if graph.edges[
        grandpa, father, LARGE
    ]['strand'][father] == '+' else operator.gt

    for adj, meta in graph.adj[father].items():
        for key in meta.keys():
            if key == PEAK and validater(adj, father):
                # larger peak first
                yield meta[key]['weight'], abs(father.start - adj.end), adj


def find_large_child(graph, grandpa, father):
    """
    find next node which is connected by a large connection

    Args:
        graph (networkx.MultiGraph): graph with large insert enriched regions as nodes
        grandpa (tuple): grandpa of the finding child
        father (tuple): father of the finding child

    Yields:
        tuple: children
    """
    direction = '+' if grandpa < father else '-'
    for adj, meta in graph.adj[father].items():
        for key in meta.keys():
            if key == LARGE and graph.edges[father, adj, LARGE]['strand'][father] == direction:
                # near connection first
                yield meta[key]['weight'], -abs(father.start - adj.end), adj


def find_children(graph, grandpa, key, father):
    """
    find next node

    Args:
        graph (networkx.MultiGraph): graph with large insert enriched regions as nodes
        grandpa (tuple): grandpa of the finding child
        key (tuple): connection type between grandpa and father
        father (tuple): father of the finding child

    Returns:
        iterable: children of given father
    """
    finder = find_peak_child if key == LARGE else find_large_child
    return (child for *_, child in sorted(finder(graph, grandpa, father), reverse=True))


def sorted_peaks(graph):
    """
    Sorted peaks in the graph

    Args:
        graph (networkx.MultiGraph): graph with large insert enriched regions as nodes

    Returns:
        iterable: peaks sorted by weight and length
    """
    peaks = []
    for (node1, node2, key), meta in graph.edges.items():
        if key == PEAK:
            node1, node2 = sorted((node1, node2))
            peaks.append(
                (meta['weight'], node2.end - node1.start, (node1, node2))
            )
    peaks.sort(reverse=True)
    return (peak for *_, peak in peaks)


def find_circle(graph):
    """
    find all possible circle in the graph by NP.

    Args:
        graph (networkx.MultiGraph): graph with large insert enriched regions as nodes

    Yields:
        tuple: unique circle
    """
    peaks = sorted_peaks(graph)

    visited = set()

    while True:
        try:
            peak = next(peaks)
        except StopIteration:
            break

        # skip visited peak as a source
        if peak in visited:
            continue
        else:
            visited.add(peak)

        circle_stack = list(peak)
        link_stack = [PEAK, ]
        stack_point = {node: i for i, node in enumerate(circle_stack)}

        stack = [
            find_children(graph, circle_stack[-2], PEAK, circle_stack[-1])
        ]

        while stack:
            children = stack[-1]

            try:
                child = next(children)
            except StopIteration:
                pop = circle_stack.pop()
                stack_point.pop(pop)
                link_stack.pop()
                stack.pop()
                continue

            if child not in stack_point:  # new peace, add to stack
                stack_point[child] = len(circle_stack)
                circle_stack.append(child)
                link_stack.append(NEXT_CONNECTION_KEY[link_stack[-1]])
                stack.append(
                    find_children(
                        graph, circle_stack[-2], link_stack[-1], child
                    )
                )
                # remember visited peaks to skip in future
                if NEXT_CONNECTION_KEY[link_stack[-1]] == PEAK:
                    peak = tuple(sorted((circle_stack[-1], child)))
                    visited.add(peak)
            elif NEXT_CONNECTION_KEY[link_stack[-1]] == link_stack[stack_point[child]]:
                continue
            else:
                if NEXT_CONNECTION_KEY[link_stack[-1]] == LARGE:
                    validator = operator.gt if graph.edges[
                        circle_stack[-1], child, LARGE
                    ]['strand'][child] == '+' else operator.lt
                    if not validator(child, circle_stack[stack_point[child] + 1]):
                        continue
                else:
                    direction = '+' if circle_stack[-1] < child else '-'
                    if graph.edges[child, circle_stack[stack_point[child] + 1], LARGE]['strand'][child] != direction:
                        continue
                yield circle_stack[stack_point[child]:], link_stack[stack_point[child]:] + [NEXT_CONNECTION_KEY[link_stack[-1]]]


def ecDNA(graph):
    """
    find ecDNA from the graph

    Args:
        graph (networkx.MultiGraph): graph with large insert enriched regions as nodes

    Yields:
        tuple: segments in a ecDNA
    """
    visited = set()
    for circle, links in find_circle(graph):
        # skip visited circle
        sc = sorted_circle(circle)
        if sc in visited:
            continue
        visited.add(sc)
        # rotate and extract segments form the circle
        segments = []
        circle, links = deque(circle), deque(links)
        begin = circle[0]
        while True:
            if links[0] == PEAK:
                segments.append((circle[0], circle[1]))
            circle.rotate(-1)
            links.rotate(-1)
            if circle[0] == begin:
                break
        yield segments


def fetch_large_pos(bam, chrom, start, end, direction, include=1, exclude=1036, mapq=10, ratio=0.05, min_insert_size=1500):
    """fetch reads positions which are within the large insert size in a specified region

    Args:
        bam (pysam.AlignmentFile): opened bam file
        chrom (string): contig name
        start (int): start
        end (int): end
        direction (string): break direction, choices from '+' or '-'
        include (int, optional): only include reads within these flags. Defaults to 1.
        exclude (int, optional): exclude reads within these flags. Defaults to 1036.
        mapq (int, optional): minimum mapq for accepted reads. Defaults to 10.
        ratio (float, optional): maximum ratio of mismatch. Defaults to 10.
        min_insert_size (int, optional): minimum length for insert size that reads pair's insert size are large. Defaults to 1500.

    Returns:
        [np.ndarray]: large insert size reads end positions
    """
    large_pos = []
    for read in bam.fetch(chrom, start, end):
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
        if read.reference_name == read.next_reference_name and abs(read.template_length) <= min_insert_size:
            continue
        if read.is_reverse == DIRECTION_FILTER[direction]:
            continue
        large_pos.append(getattr(read, POS_ATTR[direction]))
    return np.array(large_pos)


def fetch_split_pos(bam, chrom, start, end, direction, include=1, exclude=1036, mapq=10, ratio=0.05):
    """fetch reads positions which are split reads in specified region.

    Args:
        bam (pysam.AlignmentFile): opened bam file
        chrom (string): contig name
        start (int): start
        end (int ): end
        direction (string): break direction
        include (int, optional): only include reads within these flags. Defaults to 1.
        exclude (int, optional): exclude reads within these flags. Defaults to 1036.
        mapq (int, optional): minimum mapq for accepted reads. Defaults to 10.
        ratio (float, optional): maximum ratio of mismatch. Defaults to 10.

    Returns:
        [np.ndarray]: large insert size reads end positions
    """
    split_pos = []
    for read in bam.fetch(chrom, start, end):
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
        if read.cigartuples[DIRECTION_INDEX[direction]][0] not in {4, 5}:
            continue
        split_pos.append(getattr(read, POS_ATTR[direction]))
    return np.array(split_pos)


def large_log_pdf(sample, hypo):
    """calculate log pdf of given hypothesis using german tank model

    Args:
        sample ([np.ndarray]): observe samples
        hypo ([np.ndarray]): hypothesis

    Returns:
        [np.ndarray]: log pdf for given hypothesis
    """
    # calculate likelihood of all hypothesis
    likelihood = np.tile(1 / hypo, (len(hypo), 1))
    # add really small value to avoid log(0)
    likelihood[np.tril_indices_from(likelihood, k=-1)] = 0
    # cumulative product of likelihood from observe sample
    s = np.subtract.outer(sample, hypo)
    s[np.all(s < 0, axis=1), 0] = 0
    _, likelihood_index = np.where(s == 0)
    likelihood = likelihood[likelihood_index]
    pmf = np.ones_like(hypo, dtype=np.float64)
    for l in likelihood:
        pmf *= l
        pmf /= np.sum(pmf)
    return np.log(pmf + SMALL_VALUE)


def split_log_pdf(sample, hypo, times=1000):
    """calculate pdf of given hypo using bootstrap model

    Args:
        sample (np.ndarray): observe sample
        hypo (np.ndarray): hypothesis
        times (int, optional): bootstrap times. Defaults to 1000.

    Returns:
        np.ndarray: log pdf for given hypothesis
    """
    pdf = np.full_like(hypo, -10, dtype=np.float64)
    if len(sample) == 0:
        return pdf
    bootstrap_sample = np.random.choice(sample, len(sample) * times)
    bootstrap_median = np.median(
        bootstrap_sample.reshape(len(sample), times), axis=0
    )
    pos, count = np.unique(bootstrap_median, return_counts=True)
    log_freq = np.log(count / np.sum(count))
    pos_index, x_index = np.where(np.equal.outer(pos, hypo))
    pdf[x_index] = log_freq[pos_index]
    return pdf


@lru_cache(maxsize=1024)
def estimate_breakpoint(bam, chrom, start, end, direction, include=1, exclude=1036, mapq=10, ratio=0.05, min_insert_size=1500, times=1000):
    """estimate break point for the given region and direction.

    Args:
        bam (pysam.AlignmentFile): opened bam file
        chrom (string): contig name
        start (int): start
        end (int ): end
        direction (string): break direction
        include (int, optional): only include reads within these flags. Defaults to 1.
        exclude (int, optional): exclude reads within these flags. Defaults to 1036.
        mapq (int, optional): minimum mapq for accepted reads. Defaults to 10.
        ratio (float, optional): maximum ratio of mismatch. Defaults to 10.
        min_insert_size (int, optional): minimum length for insert size that reads pair's insert size are large. Defaults to 1500.
        times (int, optional): bootstrap times. Defaults to 1000.

    Returns:
        [type]: [description]
    """
    with pysam.AlignmentFile(bam, 'rb') as f:
        split_pos = fetch_split_pos(
            f, chrom, start, end, direction, include, exclude, mapq, ratio
        )
        large_pos = fetch_large_pos(
            f, chrom, start, end, direction, include, exclude, mapq, ratio, min_insert_size
        )
    offset = OFFSET_FUNC[direction](large_pos)
    large_sample = np.unique(np.abs(large_pos - offset)) + 1
    split_sample = np.abs(
        split_pos[SPLIT_POS_FILTER[direction](split_pos, offset)] - offset
    ) + 1
    hypo = np.arange(0, min_insert_size) + 1
    abs_hypo = HYPO_TRANS[direction](offset, hypo) - 1
    log_pdf = large_log_pdf(large_sample, hypo) + \
        split_log_pdf(split_sample, hypo, times)
    mle = abs_hypo[np.argmax(log_pdf)]
    cdf = np.cumsum(np.exp(log_pdf)) / np.sum(np.exp(log_pdf))
    cl, cr = np.sort(
        abs_hypo[np.sum(np.stack([cdf < 0.025, cdf < 0.975]), axis=1)]
    )
    return cl, mle, cr


def estimate_segment(region1, region2, bam, include=1, exclude=1036, mapq=10, ratio=0.05, min_insert_size=1500, times=1000):
    if region1 < region2:
        direction1, direction2 = '-', '+'
        strand = '+'
    else:
        direction1, direction2 = '+', '-'
        strand = '-'
    cl1, mle1, cr1 = estimate_breakpoint(
        bam, *region1, direction1, include, exclude, mapq, ratio, min_insert_size, times
    )
    cl2, mle2, cr2 = estimate_breakpoint(
        bam, *region2, direction2, include, exclude, mapq, ratio, min_insert_size, times
    )
    return region1[0], mle1, mle2, strand, f'{cl1}-{cr1}', f'{cl2}-{cr2}'


def run(peaks, bam, out, min_depth, include=1, exclude=1036, mapq=10, ratio=0.05, min_insert_size=1500, times=1000, limit=1000):
    """
    process to find all possible eccDNA from bam and peak files

    Args:
        peaks (str): path to the peaks file
        bam (str): path to the bam file
        out (str): path to the output file
        min_depth (int): the minimun depth of connection
        include (int, optional): only include reads within these flags. Defaults to 1.
        exclude (int, optional): exclude reads within these flags. Defaults to 1036.
        mapq (int, optional): minimum mapq for accepted reads. Defaults to 10.
        ratio (float, optional): maximum ratio of mismatch. Defaults to 10.
        min_insert_size (int, optional): minimum length for insert size that reads pair's insert size are large. Defaults to 1500.
        times (int, optional): bootstrap times. Defaults to 1000.
        limit (int, optional): limit of output
    """
    # build a graph that within large insert enriched regions
    graph = peaks2graph(peaks)
    # fetch paired reads which aligned in two large insert enriched regions
    # which means these two regions is a continue sequence in the ecDNA
    with pysam.AlignmentFile(bam, 'rb') as f:
        graph = add_large_connections(
            graph, f, min_depth, include, exclude, mapq, ratio, min_insert_size
        )
    # simplify the graph and remove nodes that cant form a circle
    nodes = [node for node, degree in graph.degree if degree > 0]
    graph = graph.subgraph(nodes).copy()
    # add peak connections
    graph = add_peak_connections(graph)
    # filter out nodes cant form a circle
    nodes = [
        node for node in graph.nodes if len(set(connect_types(graph, node))) > 1
    ]
    graph = graph.subgraph(nodes).copy()
    # find all validate ecDNA circles
    circles = ecDNA(graph)
    # output
    with open(out, 'w') as f:
        for n, circle in enumerate(circles, start=1):
            if n > 0 and n > limit:
                break
            for p, (region1, region2) in enumerate(circle, start=1):
                chrom, start, end, strand, start_ci, end_ci = estimate_segment(
                    region1, region2, bam, include, exclude, mapq, ratio, min_insert_size, times
                )
                print(
                    f'{chrom}\t{start}\t{end}\tecDNA_{n}_{p}\t.\t{strand}\t{start_ci}\t{end_ci}', file=f
                )


def main():
    parser = ArgumentParser(
        description='Find ecDNA from bam file and peaks'
    )
    parser.add_argument(
        '-d', dest='depth', help='minimum depth of connections', required=True, type=int
    )
    parser.add_argument(
        '-l', dest='length', help='minimum length of insert size', default=1500, type=int
    )
    parser.add_argument(
        '-f', dest='include', help='only include reads with all  of the FLAGs in INT present',
        default=1, type=int
    )
    parser.add_argument(
        '-F', dest='exclude', help='only include reads with none of the FLAGS in INT present',
        default=1036, type=int
    )
    parser.add_argument(
        '-q', dest='mapq', help='only include reads with mapping quality >= INT', default=10, type=int
    )
    parser.add_argument(
        '-r', dest='ratio', help='only include reads with mismatch ratio <= FLOAT', default=0.05, type=float
    )
    parser.add_argument(
        '-t', dest='times', help='bootstrap times', default=1000, type=int
    )
    parser.add_argument(
        '-m', dest='limit', help='maximum output circle.'
        ' Warning: For some heavily rearranged samples,'
        ' one can form massive possible connections between segments,'
        ' this can avoid massive output and long time running.'
        ' set to 0 to get all output (which may not be necessary).',
        default=1000, type=int
    )

    parser.add_argument(
        '<peaks>', help='large insert size reads enrich regions'
    )
    parser.add_argument('<bam>', help='input bam file')
    parser.add_argument('<bed>', help='output bed file')
    args = vars(parser.parse_args())
    run(
        args['<peaks>'], args['<bam>'], args['<bed>'], args['depth'],
        args['include'], args['exclude'], args['mapq'], args['ratio'],
        args['length'], args['times'], args['limit']
    )


if __name__ == '__main__':
    main()
