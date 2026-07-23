#!/usr/bin/env python3

"""
Command-line interface to spacetrees' core inference functions,
estimate_dispersal and locate_ancestors.

Expects the per-locus files produced by the Snakefile rules
loci_positions -> sample_trees -> extract_times -> process_times
(the .stss, .stss_logdet, _stss_inv.npy, .btss, and .lpcs files), plus a
sample locations file. This is an alternative to driving those two rules
through Snakemake, e.g. for scripting or interactive use.
"""

import argparse

import numpy as np

from spacetrees import (
    estimate_dispersal,
    locate_ancestors,
    _log_birth_density,
    _sds_rho_to_sigma,
)
from utils import chop_shared_times


def _load_stss_inv(path, m):
    """Load a *_stss_inv.npy file: M vectorized (upper triangle + diagonal) m x m matrices, one per sampled tree."""
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


def cmd_estimate_dispersal(args):

    locations = np.loadtxt(args.locations)
    n = locations.shape[0]
    m = n - 1  # samples remaining after mean-centering

    stss_logdet = [np.loadtxt(f) for f in args.stss_logdet]
    stss_inv = [_load_stss_inv(f, m) for f in args.stss_inv]

    important = not args.no_importance
    btss = None
    lpcs = None
    if important:
        if not args.btss or not args.lpcs:
            raise SystemExit('estimate-dispersal: --btss and --lpcs are required unless --no-importance is set')
        if len(args.btss) != len(args.stss_logdet) or len(args.lpcs) != len(args.stss_logdet):
            raise SystemExit('estimate-dispersal: --stss-logdet, --stss-inv, --btss, and --lpcs must all list the same number of loci, in the same order')
        btss = [_load_btss(f) for f in args.btss]
        lpcs = [np.loadtxt(f) for f in args.lpcs]

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

    # shared times (chopped at T, in matrix form), one per sampled tree
    stss_raw = np.loadtxt(args.stss, delimiter=',')
    if stss_raw.ndim == 1:
        stss_raw = stss_raw[None, :]
    stss = []
    for sts in stss_raw:
        sts = chop_shared_times(sts, T=args.T)
        mat = np.zeros((n, n))
        mat[np.triu_indices(n, k=1)] = sts[1:]
        mat = mat + mat.T + np.diag([sts[0]] * n)
        stss.append(mat)

    # shared times, mean-centered and inverted
    stss_inv = _load_stss_inv(args.stss_inv, m)

    # branching times and log coalescent-time probabilities, for importance weights
    btss = _load_btss(args.btss)
    lpcs = np.loadtxt(args.lpcs)

    sigma_raw = np.loadtxt(args.sigma, delimiter=',')
    phi = sigma_raw[-1]  # branching rate is always the last entry (see estimate-dispersal)
    sigma = _sds_rho_to_sigma(sigma_raw[:-1])

    lbds = np.array([_log_birth_density(bts, phi, n) for bts in btss])
    log_weights = lbds - lpcs

    samples = range(n) if args.samples is None else args.samples

    ancestor_locations = locate_ancestors(
        samples=samples,
        times=args.times,
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
        description="Estimate dispersal rates and locate genetic ancestors from preprocessed genealogies "
                    "(the output of the Snakefile's process_times rule).",
    )
    sub = parser.add_subparsers(dest='command', required=True)

    p1 = sub.add_parser('estimate-dispersal', help='estimate the maximum likelihood dispersal (and branching) rate')
    p1.add_argument('--stss-logdet', nargs='+', required=True, metavar='FILE', help='one *.stss_logdet file per locus')
    p1.add_argument('--stss-inv', nargs='+', required=True, metavar='FILE', help='one *_stss_inv.npy file per locus, same order as --stss-logdet')
    p1.add_argument('--btss', nargs='+', metavar='FILE', help='one *.btss file per locus, same order as --stss-logdet (required unless --no-importance)')
    p1.add_argument('--lpcs', nargs='+', metavar='FILE', help='one *.lpcs file per locus, same order as --stss-logdet (required unless --no-importance)')
    p1.add_argument('--locations', required=True, metavar='FILE', help='sample locations file')
    p1.add_argument('--out', required=True, metavar='FILE', help='where to write the estimated rate(s) (comma-separated: sdx[,sdy,rho][,phi])')
    p1.add_argument('--no-importance', action='store_true', help="don't importance-sample over branching times; ignores --btss/--lpcs and only estimates dispersal")
    p1.add_argument('--sigma0', type=float, nargs='+', default=None, metavar='X', help='initial guess, as sdx [sdy rho]; default: estimated automatically from the data')
    p1.add_argument('--quiet', action='store_true')
    p1.set_defaults(func=cmd_estimate_dispersal)

    p2 = sub.add_parser('locate-ancestors', help='locate the genetic ancestor of one or more samples at one or more times, at a single locus')
    p2.add_argument('--stss', required=True, metavar='FILE', help='*.stss file for the locus')
    p2.add_argument('--stss-inv', required=True, metavar='FILE', help='*_stss_inv.npy file for the locus (must match --T)')
    p2.add_argument('--btss', required=True, metavar='FILE', help='*.btss file for the locus (must match --T)')
    p2.add_argument('--lpcs', required=True, metavar='FILE', help='*.lpcs file for the locus (must match --T)')
    p2.add_argument('--locations', required=True, metavar='FILE', help='sample locations file')
    p2.add_argument('--sigma', required=True, metavar='FILE', help='dispersal (and branching) rate file, from estimate-dispersal')
    p2.add_argument('--out', required=True, metavar='FILE', help='where to write the located ancestors')
    p2.add_argument('--samples', type=int, nargs='+', default=None, metavar='I', help='0-indexed sample(s) to locate ancestors for; default: all samples')
    p2.add_argument('--times', type=float, nargs='+', required=True, metavar='T', help='time(s) in the past to locate ancestors at')
    p2.add_argument('--T', type=float, default=None, metavar='T', help='time cutoff the shared times were processed with (must match the --stss-inv/--btss/--lpcs files); default: none')
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
