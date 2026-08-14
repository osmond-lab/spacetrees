#!/usr/bin/env python3

"""
Command-line interface to spacetrees' core inference functions,
estimate_dispersal and locate_ancestors, plus the tree-processing steps
(loci-positions, extract-times, process-times) that produce their inputs.

Covers the pipeline
loci-positions -> sample_trees -> extract-times -> process-times -> estimate-dispersal/locate-ancestors
(the .loci, .stss, .ctss, .stss-logdet, _stss-inv.npy, .btss, and .lpcs files), plus a
sample locations file, starting from a Relate .mut file (sample_trees, run via
Relate's SampleBranchLengths.sh, is not wrapped here).
"""

import argparse
import contextlib

import numpy as np
from tqdm import tqdm

from spacetrees import (
    estimate_dispersal,
    locate_ancestors,
    _log_birth_density,
    _sds_rho_to_sigma,
)
from utils import chop_shared_times, center_shared_times, get_shared_times, loci_positions, log_coal_density


def _load_stss_inv(path, m):
    """Load a *_stss-inv.npy file: M vectorized (upper triangle + diagonal) m x m matrices, one per sampled tree."""
    vectors = np.load(path)
    matrices = []
    for vec in vectors:
        mat = np.zeros((m, m))
        mat[np.triu_indices(m, k=0)] = vec
        mat = mat + mat.T - np.diag(np.diag(mat))
        matrices.append(mat)
    return matrices


def _load_btss(path):
    """Load a *.btss file: one comma-separated list of branching times per sampled tree (one per line)."""
    btss = []
    with open(path, 'r') as f:
        for line in f:
            btss.append(np.fromstring(line, dtype=float, sep=','))
    return btss


def _T_label(T):
    """Format a time cutoff as the '{T}T' suffix used in process-times output filenames."""
    return (str(int(T)) if float(T).is_integer() else str(T)) + 'T'


def cmd_loci_positions(args):
    """Writes the position of the first and last mutation at each locus in a
    Relate .mut file, one locus per line, space-delimited."""

    loci_positions(args.mut, args.out)
    print(f'wrote loci positions to {args.out}')


def cmd_extract_times(args):
    """Pulls shared times between every pair of samples and coalescence times
    out of each sampled tree at a locus (a newick file produced by Relate's
    SampleBranchLengths.sh)."""

    from tsconvert import from_newick  # imported lazily: not needed by the other subcommands

    out_stss = args.out + '.stss'
    out_ctss = args.out + '.ctss'

    with open(args.newick, 'r') as f:
        with open(out_stss, 'w') as stss_out:
            with open(out_ctss, 'w') as ctss_out:

                next(f)  # skip header
                lines = f if args.quiet else tqdm(f)
                for line in lines:

                    string = line.split()[4]  # extract newick string only (Relate adds some info beforehand)
                    ts = from_newick(string, min_edge_length=1e-6)
                    tree = ts.first()

                    samples = [int(ts.node(node).metadata['name']) for node in ts.samples()]  # index of each sample in list given to relate
                    sample_order = np.argsort(samples)
                    ordered_samples = [ts.samples()[i] for i in sample_order]
                    sts = get_shared_times(tree, ordered_samples)  # shared times between all pairs, ordered as in relate
                    stss_out.write(",".join(str(i) for i in sts) + '\n')

                    cts = sorted([tree.time(i) for i in tree.nodes() if not tree.is_sample(i)])  # coalescence times, ascending
                    ctss_out.write(",".join(str(i) for i in cts) + '\n')

    print(f'wrote shared times to {out_stss} and coalescence times to {out_ctss}')


def cmd_process_times(args):
    """Chops shared times at a time cutoff T, centers and inverts them, and
    derives branching times and coalescent-time log probabilities, for every
    sampled tree at a locus."""

    T = args.T

    in_stss = args.times + '.stss'
    in_ctss = args.times + '.ctss'

    # with no cutoff, chopping is a no-op, so PREFIX.stss would just duplicate extract-times'
    # output: skip it and don't suffix the prefix (reuse --times' prefix for --out to let
    # locate-ancestors find that unchopped .stss alongside these processed files)
    out_prefix = args.out if T is None else args.out + '_' + _T_label(T)
    out_stss = None if T is None else out_prefix + '.stss'
    out_stss_logdet = out_prefix + '.stss-logdet'
    out_stss_inv = out_prefix + '_stss-inv.npy'
    out_btss = out_prefix + '.btss'
    out_lpcs = out_prefix + '.lpcs'

    epochs = np.genfromtxt(args.coal, skip_header=1, skip_footer=1)  # time each epoch starts (and the final one ends)
    Nes = 0.5 / np.genfromtxt(args.coal, skip_header=2)[2:]  # effective population size during each epoch

    sts_inv = []
    with open(in_stss, 'r') as stss_in:
        with open(in_ctss, 'r') as ctss_in:
            with contextlib.ExitStack() as stack:
                stss_out = stack.enter_context(open(out_stss, 'w')) if out_stss is not None else None
                stss_logdet_out = stack.enter_context(open(out_stss_logdet, 'w'))
                btss_out = stack.enter_context(open(out_btss, 'w'))
                lpcs_out = stack.enter_context(open(out_lpcs, 'w'))

                pairs = zip(stss_in, ctss_in)
                pairs = pairs if args.quiet else tqdm(pairs)
                for sts, cts in pairs:

                    # shared times: chop, convert to matrix, center, take determinant and inverse
                    sts = np.fromstring(sts, dtype=float, sep=',')
                    sts = chop_shared_times(sts, T=T)
                    if stss_out is not None:
                        stss_out.write(",".join(str(i) for i in sts) + '\n')  # same vectorized form as extract-times' *.stss

                    k = int((np.sqrt(1 + 8 * (len(sts) - 1)) + 1) / 2)  # number of samples
                    sts_mat = np.zeros((k, k))
                    sts_mat[np.triu_indices(k, k=1)] = sts[1:]
                    sts_mat = sts_mat + sts_mat.T + np.diag([sts[0]] * k)
                    sts = center_shared_times(sts_mat)

                    sts_logdet = np.linalg.slogdet(sts)[1]
                    stss_logdet_out.write(str(sts_logdet) + '\n')

                    sts = np.linalg.inv(sts)
                    sts_inv.append(sts[np.triu_indices(k - 1, k=0)])  # vectorized, to save as one array

                    # coalescence times: derive branching times and their log probability
                    cts = np.fromstring(cts, dtype=float, sep=',')
                    Tmax = cts[-1]  # time to most recent common ancestor
                    if T is not None and T < Tmax:
                        Tmax = T
                    bts = Tmax - np.flip(cts)  # branching times, ascending
                    bts = bts[bts > 0]
                    bts = np.append(bts, Tmax)  # append total time as last item
                    btss_out.write(",".join(str(i) for i in bts) + '\n')

                    lpc = log_coal_density(times=cts, Nes=Nes, epochs=epochs, T=Tmax)
                    lpcs_out.write(str(lpc) + '\n')

    with open(out_stss_inv, 'wb') as f:
        np.save(f, np.array(sts_inv))  # numpy array to avoid the numerical errors text output would introduce

    written = [out_stss, out_stss_logdet, out_stss_inv, out_btss, out_lpcs]
    print('wrote processed times to ' + ', '.join(p for p in written if p is not None))


def cmd_estimate_dispersal(args):

    locations = np.loadtxt(args.locations)
    n = locations.shape[0]
    m = n - 1  # samples remaining after mean-centering

    stss_logdet = [np.loadtxt(prefix + '.stss-logdet') for prefix in args.in_]
    stss_inv = [_load_stss_inv(prefix + '_stss-inv.npy', m) for prefix in args.in_]

    important = not args.no_importance
    btss = None
    lpcs = None
    if important:
        btss = [_load_btss(prefix + '.btss') for prefix in args.in_]
        lpcs = [np.loadtxt(prefix + '.lpcs') for prefix in args.in_]

    def callbackF(x):
        print(' '.join(f'{v: .6f}' for v in x))

    sigma = estimate_dispersal(
        locations=locations,
        shared_times_inverted=stss_inv,
        shared_times_logdet=stss_logdet,
        sigma0=args.sigma0,
        important=important,
        branching_times=btss,
        logpcoals=lpcs,
        quiet=args.quiet,
        callbackF=None if args.quiet else callbackF,
    )

    with open(args.out, 'w') as f:
        f.write(','.join(str(i) for i in sigma))
    print(f'wrote estimated rate(s) to {args.out}')


def cmd_locate_ancestors(args):

    locations = np.loadtxt(args.locations)
    n = locations.shape[0]
    m = n - 1  # samples remaining after mean-centering

    # shared times (already chopped at T by process-times, in matrix form), one per sampled tree
    stss_raw = np.loadtxt(args.in_ + '.stss', delimiter=',')
    if stss_raw.ndim == 1:
        stss_raw = stss_raw[None, :]
    stss = []
    for sts in stss_raw:
        mat = np.zeros((n, n))
        mat[np.triu_indices(n, k=1)] = sts[1:]
        mat = mat + mat.T + np.diag([sts[0]] * n)
        stss.append(mat)

    # shared times, mean-centered and inverted
    stss_inv = _load_stss_inv(args.in_ + '_stss-inv.npy', m)

    # branching times and log coalescent-time probabilities, for importance weights
    btss = _load_btss(args.in_ + '.btss')
    lpcs = np.loadtxt(args.in_ + '.lpcs')

    sigma_raw = np.loadtxt(args.sigma, delimiter=',')
    phi = sigma_raw[-1]  # branching rate is always the last entry (see estimate-dispersal)
    sigma = _sds_rho_to_sigma(sigma_raw[:-1])

    lbds = np.array([_log_birth_density(bts, phi, n) for bts in btss])
    log_weights = lbds - lpcs

    samples = range(n) if args.samples is None else args.samples

    ancestor_locations = locate_ancestors(
        samples=samples,
        times=args.ancestor_times,
        shared_times_chopped=stss,
        shared_times_chopped_centered_inverted=stss_inv,
        locations=locations,
        sigma=sigma,
        log_weights=log_weights,
        BLUP=args.blup,
        BLUP_var=args.blup_var,
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

    p0a = sub.add_parser('extract-times', help='extract shared and coalescence times from sampled trees at a locus (a newick file from Relate)')
    p0a.add_argument('--newick', required=True, metavar='FILE', help='newick file of sampled trees at a locus, as produced by SampleBranchLengths.sh')
    p0a.add_argument('--out', required=True, metavar='PREFIX', help='prefix for the output files: shared times are written to PREFIX.stss, coalescence times to PREFIX.ctss')
    p0a.add_argument('--quiet', action='store_true')
    p0a.set_defaults(func=cmd_extract_times)

    p0b = sub.add_parser('process-times', help='chop/center/invert shared times and derive branching times and coalescent-time log probabilities, at a given time cutoff T')
    p0b.add_argument('--times', required=True, metavar='PREFIX', help='prefix for the input times files, from extract-times: reads PREFIX.stss and PREFIX.ctss')
    p0b.add_argument('--coal', required=True, metavar='FILE', help="Relate's *.coal file (effective population size through time)")
    p0b.add_argument('--T', type=float, default=None, metavar='T', help='time cutoff, ignoring history beyond this time; default: none')
    p0b.add_argument('--out', required=True, metavar='PREFIX', help="base prefix for the output files. With --T, the cutoff is appended (e.g. PREFIX_10000T.stss-logdet, PREFIX_10000T_stss-inv.npy, PREFIX_10000T.btss, PREFIX_10000T.lpcs, PREFIX_10000T.stss). With no --T, chopping is a no-op, so no .stss is written and files go straight to PREFIX.stss-logdet/_stss-inv.npy/.btss/.lpcs — give --out the same prefix you gave --times so PREFIX.stss (from extract-times) sits alongside them")
    p0b.add_argument('--quiet', action='store_true')
    p0b.set_defaults(func=cmd_process_times)

    p1 = sub.add_parser('estimate-dispersal', help='estimate the maximum likelihood dispersal (and branching) rate')
    p1.add_argument('--in', dest='in_', nargs='+', required=True, metavar='PREFIX', help='one prefix per locus, from process-times: reads PREFIX.stss-logdet, PREFIX_stss-inv.npy, PREFIX.btss, and PREFIX.lpcs (the last two unused with --no-importance)')
    p1.add_argument('--locations', required=True, metavar='FILE', help='sample locations file')
    p1.add_argument('--out', required=True, metavar='FILE', help='where to write the estimated rate(s) (comma-separated: sdx[,sdy,rho][,phi])')
    p1.add_argument('--no-importance', action='store_true', help="don't importance-sample over branching times; ignores the .btss/.lpcs files and only estimates dispersal")
    p1.add_argument('--sigma0', type=float, nargs='+', default=None, metavar='X', help='initial guess, as sdx [sdy rho]; default: estimated automatically from the data')
    p1.add_argument('--quiet', action='store_true')
    p1.set_defaults(func=cmd_estimate_dispersal)

    p2 = sub.add_parser('locate-ancestors', help='locate the genetic ancestor of one or more samples at one or more times, at a single locus')
    p2.add_argument('--in', dest='in_', required=True, metavar='PREFIX', help='prefix from process-times: reads PREFIX.stss (already chopped at T), PREFIX_stss-inv.npy, PREFIX.btss, and PREFIX.lpcs')
    p2.add_argument('--locations', required=True, metavar='FILE', help='sample locations file')
    p2.add_argument('--sigma', required=True, metavar='FILE', help='dispersal (and branching) rate file, from estimate-dispersal')
    p2.add_argument('--out', required=True, metavar='FILE', help='where to write the located ancestors')
    p2.add_argument('--samples', type=int, nargs='+', default=None, metavar='I', help='0-indexed sample(s) to locate ancestors for; default: all samples')
    p2.add_argument('--ancestor_times', type=float, nargs='+', required=True, metavar='T', help='time(s) in the past to locate ancestors at')
    p2.add_argument('--blup', action='store_true', help='use the (faster) best linear unbiased predictor instead of the full likelihood surface')
    p2.add_argument('--blup-var', action='store_true', help='also report the variance of the BLUP estimate (only used with --blup)')
    p2.add_argument('--quiet', action='store_true')
    p2.set_defaults(func=cmd_locate_ancestors)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == '__main__':
    main()
