#!/usr/bin/env python3

"""
Command-line interface to spacetrees' core inference functions,
estimate_dispersal and locate_ancestors, plus the tree-processing steps
(loci-positions, process-times) that produce their inputs.

Covers the pipeline
loci-positions -> sample_trees -> process-times -> estimate-dispersal/locate-ancestors
(the .loci and .times.pkl files), plus a sample locations file, starting from
a Relate .mut file (sample_trees, run via Relate's SampleBranchLengths.sh, is
not wrapped here).
"""

import argparse
import pickle

import numpy as np
from tqdm import tqdm

from spacetrees import estimate_dispersal, locate_ancestors
from utils import split_shared_times, center_shared_times, get_shared_times, loci_positions, log_coal_density


def _T_label(T):
    """Format a time cutoff as the '{T}T' suffix used in process-times output filenames."""
    return (str(int(T)) if float(T).is_integer() else str(T)) + 'T'


def cmd_loci_positions(args):
    """Writes the position of the first and last mutation at each locus in a
    Relate .mut file, one locus per line, space-delimited."""

    loci_positions(args.mut, args.out)
    print(f'wrote loci positions to {args.out}')


def cmd_process_times(args):
    """Extracts shared and coalescence times from each sampled tree at a locus
    (a newick file produced by Relate's SampleBranchLengths.sh), chops the
    shared times at a time cutoff T and splits them into isolated subtrees,
    centers and inverts each, and derives branching times and coalescent-time
    log probabilities."""

    from tsconvert import from_newick  # imported lazily: not needed by the other subcommands

    T = args.T

    out_prefix = args.out if T is None else args.out + '_' + _T_label(T)
    out_trees = out_prefix + '.times.pkl'  # with no cutoff, one trivial subtree (all samples) per tree

    epochs = np.genfromtxt(args.coal, skip_header=1, skip_footer=1)  # time each epoch starts (and the final one ends)
    Nes = 0.5 / np.genfromtxt(args.coal, skip_header=2)[2:]  # effective population size during each epoch

    trees_data = []
    with open(args.newick, 'r') as newick_in:

        next(newick_in)  # skip header
        lines = newick_in if args.quiet else tqdm(newick_in)
        for line in lines:

            string = line.split()[4]  # extract newick string only (Relate adds some info beforehand)
            ts = from_newick(string, min_edge_length=1e-6)
            tree = ts.first()

            samples = [int(ts.node(node).metadata['name']) for node in ts.samples()]  # index of each sample in list given to relate
            sample_order = np.argsort(samples)
            ordered_samples = [ts.samples()[i] for i in sample_order]
            sts_raw = np.array(get_shared_times(tree, ordered_samples))  # shared times between all pairs, ordered as in relate

            cts = np.array(sorted(tree.time(i) for i in tree.nodes() if not tree.is_sample(i)))  # coalescence times, ascending

            # isolated subtrees: groups of samples still sharing time since T, split
            # apart from each other (their true coalescence predates the cutoff);
            # with no cutoff this is just the one subtree of all samples
            subtree_records = []
            for sample_ids, sub_mat in split_shared_times(sts_raw, T=T):
                sub_centered = center_shared_times(sub_mat)
                subtree_records.append({
                    'sample_ids': sample_ids,
                    'shared_times': sub_mat,
                    'shared_times_logdet': np.linalg.slogdet(sub_centered)[1],
                    'shared_times_inv': np.linalg.inv(sub_centered) if sub_centered.size else sub_centered,
                })

            # coalescence times: derive branching times and their log probability
            Tmax = cts[-1]  # time to most recent common ancestor
            if T is not None and T < Tmax:
                Tmax = T
            bts = Tmax - np.flip(cts)  # branching times, ascending
            bts = bts[bts > 0]
            bts = np.append(bts, Tmax)  # append total time as last item

            lpc = log_coal_density(times=cts, Nes=Nes, epochs=epochs, T=Tmax)

            # branching times and coalescent-time log probability are properties of the
            # whole tree, not of individual subtrees, so they sit alongside 'subtrees'
            trees_data.append({
                'subtrees': subtree_records,
                'branching_times': bts,
                'logpcoal': lpc,
            })

    with open(out_trees, 'wb') as f:
        pickle.dump(trees_data, f)  # list of trees, each {subtrees, branching_times, logpcoal}

    print(f'wrote processed times to {out_trees}')


def cmd_estimate_dispersal(args):

    locations = np.loadtxt(args.locations)

    processed_times = []
    for prefix in args.in_:
        with open(prefix + '.times.pkl', 'rb') as f:
            processed_times.append(pickle.load(f))

    important = not args.no_importance

    sigma = estimate_dispersal(
        locations=locations,
        processed_times=processed_times,
        important=important,
        quiet=args.quiet,
    )

    with open(args.out, 'w') as f:
        f.write(','.join(str(i) for i in sigma))
    print(f'wrote estimated rate(s) to {args.out}')


def cmd_locate_ancestors(args):

    locations = np.loadtxt(args.locations)
    n = locations.shape[0]

    with open(args.in_ + '.times.pkl', 'rb') as f:
        trees = pickle.load(f)

    sigma = None if args.sigma is None else np.loadtxt(args.sigma, delimiter=',')

    samples = range(n) if args.samples is None else args.samples

    ancestor_locations = locate_ancestors(
        samples=samples,
        ancestor_times=args.ancestor_times,
        processed_times=trees,
        sample_locations=locations,
        sigma=sigma,
        BLUP=args.blup,
        quiet=args.quiet,
    )

    with open(args.out, 'w') as f:
        for anc_loc in ancestor_locations:
            f.write(','.join([str(int(anc_loc[0]))] + [str(i) for i in anc_loc[1:]]) + '\n')
    print(f'wrote ancestor locations to {args.out}')


def build_parser():

    parser = argparse.ArgumentParser(
        prog='spacetrees',
        description="Process genealogies and estimate dispersal rates / locate genetic ancestors.",
    )
    sub = parser.add_subparsers(dest='command', required=True)

    p0 = sub.add_parser('loci-positions', help='get the position of the first and last mutation at each locus in a Relate .mut file')
    p0.add_argument('--mut', required=True, metavar='FILE', help="Relate's *.mut file")
    p0.add_argument('--out', required=True, metavar='FILE', help='where to write the locus positions (*.loci)')
    p0.set_defaults(func=cmd_loci_positions)

    p0b = sub.add_parser('process-times', help='extract shared/coalescence times from sampled trees at a locus, chop the shared times into isolated subtrees, center/invert each, and derive branching times and coalescent-time log probabilities, at a given time cutoff T')
    p0b.add_argument('--newick', required=True, metavar='FILE', help='newick file of sampled trees at a locus, as produced by SampleBranchLengths.sh')
    p0b.add_argument('--coal', required=True, metavar='FILE', help="Relate's *.coal file (effective population size through time)")
    p0b.add_argument('--T', type=float, default=None, metavar='T', help='time cutoff, ignoring history beyond this time; default: none')
    p0b.add_argument('--out', required=True, metavar='PREFIX', help="base prefix for the output file. With --T, the cutoff is appended (e.g. PREFIX_10000T.times.pkl). PREFIX.times.pkl holds, per sampled tree, its isolated subtrees since T (sample_ids/shared_times/shared_times_logdet/shared_times_inv — one trivial all-samples subtree with no --T) plus that tree's branching_times and logpcoal")
    p0b.add_argument('--quiet', action='store_true')
    p0b.set_defaults(func=cmd_process_times)

    p1 = sub.add_parser('estimate-dispersal', help='estimate the maximum likelihood dispersal (and branching) rate')
    p1.add_argument('--in', dest='in_', nargs='+', required=True, metavar='PREFIX', help='one prefix per locus, from process-times: reads PREFIX.times.pkl')
    p1.add_argument('--locations', required=True, metavar='FILE', help='sample locations file')
    p1.add_argument('--out', required=True, metavar='FILE', help='where to write the estimated rate(s) (comma-separated: sdx[,sdy,rho][,phi])')
    p1.add_argument('--no-importance', action='store_true', help="don't importance-sample over branching times; ignores each tree's branching_times/logpcoal and only estimates dispersal")
    p1.add_argument('--quiet', action='store_true')
    p1.set_defaults(func=cmd_estimate_dispersal)

    p2 = sub.add_parser('locate-ancestors', help='locate the genetic ancestor of one or more samples at one or more times, at a single locus')
    p2.add_argument('--in', dest='in_', required=True, metavar='PREFIX', help='prefix from process-times: reads PREFIX.times.pkl')
    p2.add_argument('--locations', required=True, metavar='FILE', help='sample locations file')
    p2.add_argument('--sigma', default=None, metavar='FILE', help='dispersal (and branching) rate file, from estimate-dispersal; required unless --blup (without it, trees cannot be importance-weighted by branching rate, so they are weighted equally instead)')
    p2.add_argument('--out', required=True, metavar='FILE', help='where to write the located ancestors')
    p2.add_argument('--samples', type=int, nargs='+', default=None, metavar='I', help='0-indexed sample(s) to locate ancestors for; default: all samples')
    p2.add_argument('--ancestor_times', type=float, nargs='+', required=True, metavar='T', help='time(s) in the past to locate ancestors at')
    p2.add_argument('--blup', action='store_true', help='use the (faster) best linear unbiased predictor instead of the full likelihood surface; also reports its variance as a trailing column')
    p2.add_argument('--quiet', action='store_true')
    p2.set_defaults(func=cmd_locate_ancestors)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == '__main__':
    main()
