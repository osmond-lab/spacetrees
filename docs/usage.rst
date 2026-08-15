Usage
=====

Running the pipeline
---------------------

``spacetrees`` is driven through ``cli.py``, a command-line interface with four subcommands
covering the whole pipeline from a Relate tree sequence to dispersal rates and located ancestors: ``loci-positions``,
``process-times``, ``estimate-dispersal``, and ``locate-ancestors``. With the
environment from :doc:`installation` activated, run any of them as::

    python cli.py <subcommand> --help

to see its full list of options.

Here are the steps to run spacetrees on the test data.

1. Get the position of the first and last mutation at each locus from a Relate ``.mut`` file::

    python cli.py loci-positions \
        --mut data/test_chr1.mut \
        --out data/test_chr1.loci

2. Sample trees at a locus using `Relate <https://myersgroup.github.io/relate/>`_'s ``SampleBranchLengths.sh`` (not wrapped by ``cli.py`` —
   run it directly). It takes the first/last base-pair position of the locus, which you can read
   from the ``.loci`` file above (e.g. locus 1 is the first line)::

    start=$(awk 'NR==1 {print $1}' data/test_chr1.loci)
    stop=$(awk 'NR==1 {print $2}' data/test_chr1.loci)

    relate/scripts/SampleBranchLengths/SampleBranchLengths.sh \
        -i data/test_chr1 \
        --dist data/test_chr1.dist \
        --coal data/test.coal \
        -o data/test_chr1_locus1 \
        -m 1e-8 \
        --format n \
        --num_samples 10 \
        --first_bp $start \
        --last_bp $stop \
        --seed 1

This writes ``data/test_chr1_locus1.newick`` (plus per-locus ``.anc``/``.mut``/``.dist``/``.sites``
files), sampling 10 trees (``--num_samples``) at a mutation rate of ``1e-8`` (``-m``).

3. Extract shared and coalescence times from the sampled trees at a locus, chop the shared times at
a given time cutoff ``T`` (splitting each tree into isolated subtrees — groups of samples still
sharing time since ``T``; their true coalescence predates the cutoff, so they're treated as
independent), center and invert each subtree's shared times, and derive branching times and
coalescent-time log probabilities. Omit ``--T`` for no cutoff (one trivial subtree containing all
samples). ``--out`` is a base prefix — the cutoff is appended automatically::

    python cli.py process-times \
        --newick data/test_chr1_locus1.newick \
        --coal data/test.coal \
        --T 10000 \
        --out data/test_chr1_locus1

This writes ``data/test_chr1_locus1_10000T.times.pkl`` — a pickled list, one entry per sampled
tree, each holding:

- ``subtrees``: the tree's isolated subtrees since ``T``, each with ``sample_ids``,
  ``shared_times``, ``shared_times_logdet``, and ``shared_times_inv``
- ``branching_times``
- ``logpcoal``

Omitting ``--T`` skips the cutoff and writes ``data/test_chr1_locus1.times.pkl`` with no suffix.

4. Estimate a dispersal rate from one or more loci::

    python cli.py estimate-dispersal \
        --in data/test_chr1_locus1_10000T data/test_chr1_locus2_10000T \
        --locations data/test.locations \
        --out data/test.sigma

Each ``--in`` entry is a full prefix from ``process-times`` (including its ``T`` suffix): this
reads ``data/test_chr1_locus1_10000T.times.pkl`` and ``data/test_chr1_locus2_10000T.times.pkl``.

``--out`` is written as a single comma-separated line. For 2D locations (as in the test data) this is ``sdx,sdy,rho,phi``: the dispersal standard
deviations along each axis, their correlation ``rho``, and the branching rate ``phi`` (for 1D
locations it's just ``sdx,phi``). ``phi`` is only present when importance sampling is used (the
default) — with ``--no-importance`` it's omitted, so that file can't be used by
``locate-ancestors`` below.

5. Locate ancestors at a locus::

    python cli.py locate-ancestors \
        --in data/test_chr1_locus1_10000T \
        --locations data/test.locations \
        --sigma data/test.sigma \
        --samples 0 1 \
        --ancestor_times 10 100 1000 \
        --out data/test_chr1_locus1.locs

``--in`` reads ``data/test_chr1_locus1_10000T.times.pkl``, the file written by ``process-times``
for that prefix (already chopped at ``T=10000``).

Add ``--blup`` to ``locate-ancestors`` for the faster best linear unbiased predictor instead of the
full likelihood surface. This does not require a dispersal estimate (omit ``--sigma``), which is the slowest part of the pipeline.

Each row of ``--out`` is ``sample,time,x,y,...`` — the 0-indexed sample (matching the row order in ``--locations``), the time in the past, and one column per spatial dimension of the estimated ancestor location at that time (two columns, ``x`` and ``y``, for the 2D test data). With
``--blup``, each row gets one extra trailing column: the variance of the BLUP estimate.
