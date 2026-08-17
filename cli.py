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

from spacetrees import estimate_dispersal, locate_ancestors
from utils import loci_positions, process_times


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

    trees_data, sample_times = process_times(args.newick, args.coal, T=args.T, quiet=args.quiet)

    out_prefix = args.out if args.T is None else args.out + '_' + _T_label(args.T)
    out_trees = out_prefix + '.times.pkl'  # with no cutoff, one trivial subtree (all samples) per tree

    with open(out_trees, 'wb') as f:
        # 'trees': list of trees, each {subtrees, branching_times, logpcoal}
        # 'sample_times': each sample's own age (0 for contemporary samples), aligned to the locations file
        pickle.dump({'trees': trees_data, 'sample_times': sample_times}, f)

    print(f'wrote processed times to {out_trees}')


def cmd_estimate_dispersal(args):

    locations = np.loadtxt(args.locations)

    processed_times = []
    sample_times = None  # assumed consistent across loci (same samples, same ages); taken from the first
    for prefix in args.in_:
        with open(prefix + '.times.pkl', 'rb') as f:
            data = pickle.load(f)
        processed_times.append(data['trees'])
        if sample_times is None:
            sample_times = data['sample_times']

    important = not args.no_importance

    sigma = estimate_dispersal(
        locations=locations,
        processed_times=processed_times,
        important=important,
        sample_times=sample_times,
        BLUP=args.blup,
        quiet=args.quiet,
    )

    with open(args.out, 'w') as f:
        f.write(','.join(str(i) for i in sigma))
    print(f'wrote estimated rate(s) to {args.out}')


def cmd_locate_ancestors(args):

    locations = np.loadtxt(args.locations)
    n = locations.shape[0]

    with open(args.in_ + '.times.pkl', 'rb') as f:
        data = pickle.load(f)
    trees = data['trees']
    sample_times = data['sample_times']

    sigma = None if args.sigma is None else np.loadtxt(args.sigma, delimiter=',')

    samples = range(n) if args.samples is None else args.samples

    ancestor_locations = locate_ancestors(
        samples=samples,
        ancestor_times=args.ancestor_times,
        processed_times=trees,
        sample_locations=locations,
        sigma=sigma,
        sample_times=sample_times,
        forget_locations=args.forget_locations,
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
    p0b.add_argument('--out', required=True, metavar='PREFIX', help="base prefix for the output file. With --T, the cutoff is appended (e.g. PREFIX_10000T.times.pkl). PREFIX.times.pkl holds {'trees': ..., 'sample_times': ...} -- 'trees' is a list, one per sampled tree, of its isolated subtrees since T (sample_ids/shared_times/shared_times_logdet/shared_times_inv — one trivial all-samples subtree with no --T) plus that tree's branching_times and logpcoal; 'sample_times' is each sample's own age (0 for contemporary samples), aligned to the locations file")
    p0b.add_argument('--quiet', action='store_true')
    p0b.set_defaults(func=cmd_process_times)

    p1 = sub.add_parser('estimate-dispersal', help='estimate the maximum likelihood dispersal (and branching) rate')
    p1.add_argument('--in', dest='in_', nargs='+', required=True, metavar='PREFIX', help='one prefix per locus, from process-times: reads PREFIX.times.pkl')
    p1.add_argument('--locations', required=True, metavar='FILE', help='sample locations file')
    p1.add_argument('--out', required=True, metavar='FILE', help='where to write the estimated rate(s) (comma-separated: sdx[,sdy,rho][,phi] -- phi omitted with --no-importance or --blup)')
    p1.add_argument('--no-importance', action='store_true', help="don't importance-sample over branching times; ignores each tree's branching_times/logpcoal and only estimates dispersal")
    p1.add_argument('--blup', action='store_true', help='quick average MLE over trees rather than numerical search for max composite likelihood; never estimates phi')
    p1.add_argument('--quiet', action='store_true')
    p1.set_defaults(func=cmd_estimate_dispersal)

    p2 = sub.add_parser('locate-ancestors', help='locate the genetic ancestor of one or more samples at one or more times, at a single locus')
    p2.add_argument('--in', dest='in_', required=True, metavar='PREFIX', help='prefix from process-times: reads PREFIX.times.pkl')
    p2.add_argument('--locations', required=True, metavar='FILE', help='sample locations file')
    p2.add_argument('--sigma', default=None, metavar='FILE', help='dispersal and branching rate file, from estimate-dispersal; required unless --blup (without it, trees cannot be importance-weighted by branching rate, so they are weighted equally instead). Must include phi (i.e. from an estimate-dispersal run without --blup or --no-importance), whether or not --blup is passed here')
    p2.add_argument('--out', required=True, metavar='FILE', help='where to write the located ancestors')
    p2.add_argument('--samples', type=int, nargs='+', default=None, metavar='I', help='0-indexed sample(s) to locate ancestors for; default: all samples')
    p2.add_argument('--ancestor_times', type=float, nargs='+', default=None, metavar='T', help='time(s) in the past to locate ancestors at; default: each sample\'s own sample time (from process-times), e.g. for --forget-locations')
    p2.add_argument('--forget-locations', dest='forget_locations', action='store_true', help="mask the given --samples' own known locations (held out together, not one at a time) and predict them from shared times with everyone else under Brownian motion -- e.g. to check the model fit, or trace a migrant lineage's origin")
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
