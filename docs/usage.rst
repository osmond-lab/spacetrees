Usage
=====

Running the pipeline
---------------------

``spacetrees`` is driven through ``cli.py``, a command-line interface with five subcommands
covering the whole pipeline from a Relate tree sequence to dispersal rates and located ancestors: ``loci-positions``,
``extract-times``, ``process-times``, ``estimate-dispersal``, and ``locate-ancestors``. With the
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

3. Extract shared and coalescence times from the sampled trees at a locus::

    python cli.py extract-times \
        --newick data/test_chr1_locus1.newick \
        --out data/test_chr1_locus1

This writes ``data/test_chr1_locus1.stss`` and ``data/test_chr1_locus1.ctss``.

4. Process those times at a given time cutoff ``T`` (chop, center, invert the shared times, and
derive branching times and coalescent-time log probabilities; omit ``--T`` for no cutoff).
``--out`` is a base prefix — the cutoff is appended automatically, so the files this writes
depend on ``--T``::

    python cli.py process-times \
        --times data/test_chr1_locus1 \
        --coal data/test.coal \
        --T 10000 \
        --out data/test_chr1_locus1

This reads ``data/test_chr1_locus1.stss`` and ``data/test_chr1_locus1.ctss`` (from ``extract-times`` above) and writes
``data/test_chr1_locus1_10000T.stss``, ``data/test_chr1_locus1_10000T.stss-logdet``, ``data/test_chr1_locus1_10000T_stss-inv.npy``,
``data/test_chr1_locus1_10000T.btss``, and ``data/test_chr1_locus1_10000T.lpcs``.

Omitting ``--T`` skips the ``.stss`` output — chopping is a no-op with no cutoff, so it would only
duplicate ``extract-times``' file — and writes the rest directly to ``--out`` with no suffix. Use
the *same* prefix as ``--times`` so the unchopped ``data/test_chr1_locus1.stss`` sits right alongside them.


5. Estimate a dispersal rate from one or more loci::

    python cli.py estimate-dispersal \
        --in data/test_chr1_locus1_10000T data/test_chr1_locus2_10000T \
        --locations data/test.locations \
        --out data/test.sigma

Each ``--in`` entry is a full prefix from ``process-times`` (including its ``T`` suffix): this
reads ``data/test_chr1_locus1_10000T.stss-logdet``, ``data/test_chr1_locus1_10000T_stss-inv.npy``, ``data/test_chr1_locus1_10000T.btss``, and
``data/test_chr1_locus1_10000T.lpcs``, and the same four files for ``data/test_chr1_locus2_10000T``.

``--out`` is written as a single comma-separated line. For 2D locations (as in the test data) this is ``sdx,sdy,rho,phi``: the dispersal standard
deviations along each axis, their correlation ``rho``, and the branching rate ``phi`` (for 1D
locations it's just ``sdx,phi``). ``phi`` is only present when importance sampling is used (the
default) — with ``--no-importance`` it's omitted, so that file can't be used directly by
``locate-ancestors`` below without also specifying ``--blup`` below.

6. Locate ancestors at a locus::

    python cli.py locate-ancestors \
        --in data/test_chr1_locus1_10000T \
        --locations data/test.locations \
        --sigma data/test.sigma \
        --samples 0 1 \
        --ancestor_times 10 100 1000 \
        --out data/test_chr1_locus1.locs

``--in`` reads all four files written by ``process-times`` for that prefix: ``data/test_chr1_locus1_10000T.stss``
(already chopped at ``T=10000``), ``data/test_chr1_locus1_10000T_stss-inv.npy``, ``data/test_chr1_locus1_10000T.btss``, and
``data/test_chr1_locus1_10000T.lpcs``.

Add ``--blup`` (and optionally ``--blup-var``) to ``locate-ancestors`` for the faster best linear
unbiased predictor instead of the full likelihood surface. This does not require a dispersal estimate (omit ``--sigma``), which is the slowest part of the pipeline.

Each row of ``--out`` is ``sample,time,x,y,...`` — the 0-indexed sample (matching the row order in ``--locations``), the time in the past, and one column per spatial dimension of the estimated ancestor location at that time (two columns, ``x`` and ``y``, for the 2D test data). With
``--blup --blup-var``, each row gets one extra trailing column: the variance of the BLUP
estimate.
