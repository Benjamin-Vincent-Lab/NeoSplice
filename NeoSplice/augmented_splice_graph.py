
"""splice_graph: Build a splice graph from a bam file of aligned reads and a reference genome.

Additional utilities are provided to enable basic manipulation of splice graphs.
"""

import argparse
from collections import defaultdict
import json
import logging
import networkx as nx
from networkx.readwrite import json_graph
import numpy as np
import os
import pickle
from pyfaidx import Fasta
import pysam
from scipy.stats import binom
import sys
from BCBio import GFF


class CigarOp:
    """Enumeration of CIGAR codes."""
    MATCH = 0
    INS = 1
    DEL = 2
    REF_SKIP = 3
    SOFT_CLIP = 4


class EdgeType:
    """Enumeration of edge types in the resulting graph."""
    EXON = 'exon'
    SPLICE = 'splice'
    DEL = 'del'
    TSS = 'tss'


class Strand:
    """Enumeration of strand types.

    Each strand corresponds to a transcription direction in the resultant graph.
    Note that (PLUS | MINUS) == BOTH. This should not be changed, as this fact
    is used in the module.
    """
    PLUS = 1
    MINUS = 2
    BOTH = 3


def adjacency_graph(data, directed=False, multigraph=True):
    """Return graph from adjacency data format.

    Parameters
    ----------
    data : dict
        Adjacency list formatted graph data

    Returns
    -------
    G : NetworkX graph
       A NetworkX graph object

    directed : bool
        If True, and direction not specified in data, return a directed graph.

    multigraph : bool
        If True, and multigraph not specified in data, return a multigraph.

    Examples
    --------
    >>> from networkx.readwrite import json_graph
    >>> G = nx.Graph([(1,2)])
    >>> data = json_graph.adjacency_data(G)
    >>> H = json_graph.adjacency_graph(data)

    See Also
    --------
    adjacency_graph, node_link_data, tree_data

    Note: This method is a fork of one in networkx.readwrite.json_graph, needed here
    to fix a defect with the current method. A NetworkX MultiDiGraph object
    includes edge keys in edge attribute dictionaries, causing a duplicate
    keyword argument error. This modification resolves that issue.

    Note: This change has not been merged into networkx master because this method has undergone significant
    changes there. The updated method may or may not have this bug.

    Source: https://github.com/networkx/networkx/blob/master/networkx/readwrite/json_graph/adjacency.py
    """
    multigraph = data.get('multigraph',multigraph)
    directed = data.get('directed',directed)
    if multigraph:
        graph = nx.MultiGraph()
    else:
        graph = nx.Graph()
    if directed:
        graph = graph.to_directed()
    graph.graph = dict(data.get('graph',[]))
    mapping=[]
    for d in data['nodes']:
        node_data = d.copy()
        node = node_data.pop('id')
        mapping.append(node)
        graph.add_node(node, attr_dict=node_data)
    for i,d in enumerate(data['adjacency']):
        source = mapping[i]
        for tdata in d:
            target_data = tdata.copy()
            target = target_data.pop('id')
            key = target_data.pop('key', None)
            if not multigraph or key is None:
                graph.add_edge(source,target,attr_dict=target_data)
            else:
                graph.add_edge(source,target,key=key, attr_dict=target_data)
    return graph


def _open_bam(bam_name):
    return pysam.AlignmentFile(bam_name, 'rb')


def find_genes(splice_graph, filter=None):
    """Find weakly connected components of the given splice graph that satisfy the
    provided filter. If no filter is provided, all compoqnents will be generated.
    """
    for i, gene in enumerate(nx.weakly_connected_component_subgraphs(splice_graph)):
        if filter is None:
            yield (i, gene)
        elif filter(gene):
            yield (i, gene)


def find_all_genes(splice_graph):
    """
    Return all genes in the splice graph
    """
    return [gene for gene in nx.weakly_connected_component_subgraphs(splice_graph)]


def _save_splice_graph_json(splice_graph, name, indent=None):
    with open('{0}.json'.format(name), 'w') as f:
        # Convert to serializable representation and write json file
        splice_graph_json = json_graph.adjacency_data(splice_graph)

        json.dump(splice_graph_json, f, indent=indent)


def save_splice_graph(splice_graph, name, format):
    """Save a splice graph in the specified format.

    If pickle is chosen, the necessary conversions will be made to support the
    legacy format.
    """
    if format == 'json':
        _save_splice_graph_json(splice_graph, name)
    else:
        logging.error('Unsupported graph output format')
        sys.exit(1)


def _load_splice_graph_json(name):
    with open(name, 'r') as f:
        # Load the json representation of the graph and convert to networkx
        splice_graph_json = json.load(f)
        return adjacency_graph(splice_graph_json)


def load_splice_graph(name):
    """Load a splice graph from the given path.

    Format is determined automatically from the file extension.
    """
    ext = os.path.splitext(name)[1]
    if ext is None or ext == '.json':
        return _load_splice_graph_json(name)
    else:
        logging.error('Unsupported graph input format')
        sys.exit(1)


def _find_TSS(gff_in_file, bam_seq):
    in_file = gff_in_file
    in_handle = open(in_file)
    limit_info = dict(gff_id=[bam_seq])
    tsss = defaultdict(list)

    for chromosome in GFF.parse(in_handle, limit_info=limit_info):
        for gene in chromosome.features:
            if 'protein_coding' in gene.qualifiers['gene_type']:
                for transcript in gene.sub_features:
                    if 'protein_coding' in transcript.qualifiers['transcript_type']:
                        for feature in transcript.sub_features:
                            feature_start = feature.location.start.position
                            feature_end = feature.location.end.position

                            if feature.type == 'start_codon':
                                tsss[transcript.id].append([feature_start, feature_end,
                                                            feature.qualifiers['transcript_id']])

    return tsss


def _find_long_variants(bam, seq, seq_len):
    # Track whether we've logged various warnings to keep the noise down
    xs_tag_warning_logged = False
    unexpected_cigar_op_warning_logged = False

    insertions = defaultdict(int)
    deletions = defaultdict(int)
    splices = defaultdict(int)
    splice_strands = defaultdict(int)
    
    for read in bam.fetch(seq, 0, seq_len):
        ref_pos = read.reference_start
        read_pos = 0

        # Get strand info, if available
        strand = None
        tags = dict(read.get_tags())
        if 'XS' in tags:
            strand_tag = tags['XS']
            if strand_tag == '+':
                strand = Strand.PLUS
            elif strand_tag == '-':
                strand = Strand.MINUS
        elif not xs_tag_warning_logged:
            xs_tag_warning_logged = True
            logging.warn('At least one splice was missing the XS tag. Both directions may be used in the graph.')

        if not read.cigartuples:
            continue

        for op, count in read.cigartuples:
            if op == CigarOp.MATCH:
                ref_pos += count
                read_pos += count
            elif op == CigarOp.INS:
                # Inserted sequence begins and ends before ref_pos
                start = ref_pos
                seq = read.query_sequence[read_pos:read_pos+count].upper()
                insertions[start, seq] += 1
                read_pos += count
            elif op == CigarOp.DEL:
                # Deletion edge extends from first deleted base to next retained base
                start = ref_pos
                end = ref_pos + count
                deletions[start, end] += 1
                ref_pos += count
            elif op == CigarOp.REF_SKIP:
                # Splice edge extends from first excised base to next retained base
                start = ref_pos
                end = ref_pos + count
                if strand is not None:
                    splice_strands[start, end] |= strand
                splices[start, end] += 1
                ref_pos += count
            elif op == CigarOp.SOFT_CLIP:
                read_pos += count
            elif not unexpected_cigar_op_warning_logged:
                unexpected_cigar_op_warning_logged = True
                logging.warn('Unexpected cigar op {0}'.format(op))

    return insertions, deletions, splices, splice_strands


def _analyze_column(column, snp_cutoff, p_error, ref_base):
    coverage = 0
    bases = defaultdict(int)
    
    for read in column.pileups:
        if not read.is_refskip and not read.is_del:
            coverage += 1
            base = read.alignment.query_sequence[read.query_position].upper()
            bases[base] += 1

    snp, snp_weight = _call_snp(bases, snp_cutoff, p_error, ref_base, coverage)
    return coverage, bases, snp, snp_weight


def _call_snp(bases, cutoff, p_error, ref_base, coverage):
    alt_base = None
    alt_count = None

    for base, count in bases.iteritems():
        if base != ref_base and (alt_count is None or count > alt_count):
            alt_base = base
            alt_count = count

    if alt_base is None:
        return None, 0

    if _p_at_least_k_errors_in_n_bases(coverage, alt_count, p_error) < cutoff:
        return alt_base, alt_count

    return None, 0


def _p_at_least_k_errors_in_n_bases(n, k, p_error):
    return binom.sf(k-1, n, p_error)


def _pileup(bam, reference, seq, seq_len, splice_graph, snp_cutoff, min_coverage, min_variants, p_error):
    start = None
    last_pos = None

    for column in bam.pileup(seq, 0, seq_len, truncate=True, max_depth=200000):
        # Analyze reads and call snps
        ref_base = reference[column.reference_pos].seq.upper()
        coverage, bases, snp, snp_weight = _analyze_column(column, snp_cutoff, p_error, ref_base)
        snps_to_add = snp is not None and snp_weight >= min_variants
        effective_coverage = coverage - snp_weight
        has_coverage = coverage >= min_coverage
        
        if start is not None:
            hit_splice = last_pos+1 in splice_graph
            skip_ahead = column.reference_pos != last_pos+1
            terminate_exon = not has_coverage or hit_splice or snps_to_add or skip_ahead

            if terminate_exon:
                # Calculate average coverage for the exon
                average_coverage = np.mean(running_coverage)

                # Get the reference sequence
                exon_seq = reference[start:last_pos+1].seq.upper()
                
                # Edge extends from first transcribed base to next base not in the exon
                splice_graph.add_edge(start, last_pos+1, EdgeType.EXON, {'weight': average_coverage, 'seq': exon_seq, 'coverage': running_coverage})

                start = None

        if has_coverage:
            if start is None:
                # If we're not in an exon or just ended one and have coverage
                running_coverage = [effective_coverage]
                start = column.reference_pos
            else:
                # Extend a transcription edge if we're in one
                running_coverage.append(effective_coverage)
        # Insert snp edge if needed. Note that the next outer loop iteration will
        # always terminate an exon if a snp edge is inserted.
        if snps_to_add:
            splice_graph.add_edge(column.reference_pos, column.reference_pos + 1, snp, {'weight': snp_weight, 'seq': snp})

        last_pos = column.reference_pos

    logging.info('Pileup finished, preparing to flush last exon. start={0}, last_pos={1}'.format(start, last_pos))

    # Handle the case where an exon was in progress on the last base examined
    if start is not None:
        # Calculate average coverage for the exon
        average_coverage = np.mean(running_coverage)

        # Get the reference sequence
        exon_seq = reference[start:last_pos+1].seq.upper()

        # Edge extends from first transcribed base to next base not in the exon
        splice_graph.add_edge(start, last_pos+1, EdgeType.EXON, {'weight': average_coverage, 'seq': exon_seq, 'coverage': running_coverage})


def _redirect_edges(splice_graph):
    # Process each connected component
    for _, gene in find_genes(splice_graph):
        edges = gene.edges(data=True, keys=True)

        # Try to figure out the direction of transcription
        strand = 0
        for edge in edges:
            # Edges are 4-tuples (u,v,key,attr_dict)
            if 'strand' in edge[3]:
                strand |= edge[3]['strand']
                if strand == Strand.BOTH:
                    break

        # Update the strand attribute on all edges in the component
        for edge in edges:
            edge[3]['strand'] = strand
            splice_graph.add_edge(edge[0], edge[1], edge[2], edge[3])


def _build_splice_graph(gff_name, bam_name, bam_seq, genome_name, genome_seq, snp_cutoff, min_coverage, min_variants, p_error):
    if genome_seq is None:
        genome_seq = bam_seq

    with _open_bam(bam_name) as bam:
        for seq_header in bam.header['SQ']:
            if seq_header['SN'] == bam_seq:
                seq_len = seq_header['LN']
                logging.info('Working on seq {0} with length {1}'.format(bam_seq, seq_len))
                break
        else:
            logging.error('Unable to find seq {0} in {1}'.format(bam_seq, bam_name))
            sys.exit(1)

        # Get the reference genome
        genome = Fasta(genome_name)
        if genome_seq not in genome:
            logging.error('Unable to find seq {0} in {1}'.format(genome_seq, genome_name))
            sys.exit(1)
        reference = genome[genome_seq]

        # Initialize the directed multi-graph
        splice_graph = nx.MultiDiGraph(name=bam_seq)
    
        # Determine multi-base variants before pileup
        logging.info('Finding insertions, deletions, and splices')
        insertions, deletions, splices, splice_strands = _find_long_variants(bam, bam_seq, seq_len)
        tsss = _find_TSS(gff_name, bam_seq)

        # add translation initiation sites to the graph
        for annotated_transcript in tsss:
            for start, end, ID in tsss[annotated_transcript]:
                splice_graph.add_edge(start, end, EdgeType.TSS, {'ID': ID})

        # Add insertion edges to the graph
        for (start, seq), count in insertions.iteritems():
            if count >= min_variants:
                splice_graph.add_edge(start, start, seq, {'weight': count, 'seq': seq})

        # Add deletion edges to the graph
        for (start, end), count in deletions.iteritems():
            if count >= min_variants:
                splice_graph.add_edge(start, end, EdgeType.DEL, {'weight': count})

        # Add splice edges to the graph
        for (start, end), count in splices.iteritems():
            strand = splice_strands[start, end]
            splice_graph.add_edge(start, end, EdgeType.SPLICE, {'weight': count, 'strand': strand})
        
        # Perform the pileup iteration, modifying splice_graph
        num_reads = bam.count(bam_seq, 0, seq_len)
        logging.info('Starting pileup iteration with {0} reads'.format(num_reads))
        _pileup(bam, reference, bam_seq, seq_len, splice_graph, snp_cutoff, min_coverage, min_variants, p_error)

        # Post-processing step
        logging.info('Starting edge redirection')
        _redirect_edges(splice_graph)

        return splice_graph


def _create_output_directory(out):
    if not os.path.isdir(out):
        logging.info('Creating output directory')
        os.mkdir(out)


def _build(args):
    splice_graph = _build_splice_graph(args.gff, args.bam, args.seq, args.genome, args.genome_seq, args.cutoff, args.min_coverage, args.min_variants, args.p_error)

    _create_output_directory(args.out)

    logging.info('Saving splice graph')
    base_output = os.path.join(args.out, args.seq)
    name = '{0}_graph'.format(base_output)
    save_splice_graph(splice_graph, name, args.format)

    logging.info('Done!')


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)-15s [%(levelname)s] %(message)s')

    parser = argparse.ArgumentParser(description='Utilities for creating and working with splice graphs.')
    subparsers = parser.add_subparsers(title='commands')

    build_parser = subparsers.add_parser('build', help='Construct a splice graph from a bam file.')
    build_parser.add_argument('--gff', required=True, help='The GFF file to use when building the graph.')
    build_parser.add_argument('--bam', required=True, help='The bam file to use when building the graph.')
    build_parser.add_argument('--seq', required=True, help='The sequence in the bam file to use.')
    build_parser.add_argument('--genome', required=True, help='The fasta file containing the reference genome.')
    build_parser.add_argument('--genome-seq', help='The sequence to use in the fasta file, if different from seq.')
    build_parser.add_argument('--p-error', type=float, default=0.0026, help='Estimate of the error probability for a base call. Defaults to 0.0026, which is calibrated for an Illumina HiSeq.')
    build_parser.add_argument('--cutoff', type=float, default=0.01, help='Significance cutoff for SNP calling. Defaults to 0.01.')
    build_parser.add_argument('--min-coverage', type=int, default=1, help='Minimum coverage to be included in the graph. Defaults to 1.')
    build_parser.add_argument('--min-variants', type=int, default=2, help='Minimum variant count to be included in the graph. Defaults to 2.')
    build_parser.add_argument('--out', default='output', help='Location to write output files. Defaults to \'output\'.')
    build_parser.add_argument('--format', default='json', help='Splice graph output format. Defaults to json.')
    build_parser.set_defaults(command=_build)

    args = parser.parse_args()
    args.command(args)


if __name__ == '__main__':
    main()
