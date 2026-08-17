Usage
=====

Running the pipeline
---------------------

``spacetrees`` is driven through ``cli.py``, a command-line interface with four subcommands
covering the whole pipeline from a `Relate <https://myersgroup.github.io/relate/>`_ geneaology
to dispersal rates and ancestor locations: ``loci-positions``,
``process-times``, ``estimate-dispersal``, and ``locate-ancestors``. With the
environment from :doc:`installation` activated, run any of them as::

    python cli.py <subcommand> --help

to see its full list of options.

Here are the steps to run spacetrees on the test data.

1. Get the position of the first and last mutation at each tree (hereafter locus) from a Relate ``.mut(.gz)`` file::

    python cli.py loci-positions \
        --mut data/test_chr1.mut \
        --out data/test_chr1.loci

2. Sample branch lengths (hereafter trees) at a locus using Relate's ``SampleBranchLengths.sh`` (not wrapped by ``cli.py`` —
   run it directly, see `Relate docs <https://myersgroup.github.io/relate/modules.html#SampleBranchLengths/>_`). 
   We feed it the first and last base-pair position of a locus, read
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
files), sampling 10 trees (``--num_samples``) assuming a mutation rate of ``1e-8`` (``-m``).

3. Extract shared and coalescence times from the sampled trees at a locus, chop the shared times at
a given time cutoff ``T`` (splitting each tree into isolated subtrees — groups of samples still
sharing time since ``T``), center and invert each subtree's shared times, and derive branching times and
coalescent-time log probabilities. Omit ``--T`` for no cutoff (one trivial subtree containing all
samples). ``--out`` is a base prefix — the cutoff is appended automatically::

    python cli.py process-times \
        --newick data/test_chr1_locus1.newick \
        --coal data/test.coal \
        --T 10000 \
        --out data/test_chr1_locus1

This writes ``data/test_chr1_locus1_10000T.times.pkl`` — a pickled dict with two entries:

- ``trees``: a list, one entry per sampled tree, each holding:

  - ``subtrees``: the tree's isolated subtrees since ``T``, each with ``sample_ids``,
    ``shared_times`` (matrix of shared times), ``shared_times_logdet`` (log determinant), and ``shared_times_inv`` (inverted shared time matrix)
  - ``branching_times``
  - ``logpcoal`` (log probability of the coalescence times given Ne in ``.coal`` file)

- ``sample_times``: each sample's own age, read straight off the genealogy

Omitting ``--T`` skips the cutoff and writes ``data/test_chr1_locus1.times.pkl`` with no suffix.

4. Estimate a dispersal rate from one or more loci::

    python cli.py estimate-dispersal \
        --in data/test_chr1_locus1_10000T data/test_chr1_locus2_10000T \
        --locations data/test.locations \
        --out data/test.sigma

Each ``--in`` entry is a full prefix from ``process-times`` (including its ``T`` suffix): this
reads ``data/test_chr1_locus1_10000T.times.pkl`` and ``data/test_chr1_locus2_10000T.times.pkl``.
``--locations`` is the location of each sample, listed in the same order as in the Relate genealogy
(i.e., as in the ``.poplabels`` file).

``--out`` is written as a single comma-separated line. For 2D locations (as in the test data)
this is ``sdx,sdy,rho,phi``: the dispersal standard
deviations along each axis, their correlation ``rho``, and the branching rate ``phi`` (for 1D
locations it's just ``sdx,phi``). ``phi`` is only present when importance sampling is used (the
default) — with ``--no-importance`` the trees are not weighted by importance and ``phi`` is not estimated.
``--blup`` (below) also never estimates ``phi``, regardless of ``--no-importance``.

Add ``--blup`` to ``estimate-dispersal`` to estimate dispersal as a weighted average of tractable maximum
likelihood over trees instead of maximizing the composite likelihood (which requires a numerical search).
Because it skips the composite-likelihood search entirely, its output never includes ``phi`` — just
``sdx[,sdy,rho]``.

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

Add ``--blup`` to ``locate-ancestors`` for the faster best linear unbiased predictor instead of the maximum
likelihood estimate. This does not require a dispersal estimate (omit ``--sigma``), which is the slowest part of the pipeline.

If you do pass ``--sigma`` (with or without ``--blup``), it must include ``phi`` — i.e. come from an
``estimate-dispersal`` run *without* ``--blup`` or ``--no-importance`` — since ``locate-ancestors`` always
uses ``phi`` to importance-weight trees by branching rate whenever a sigma is supplied. A sigma missing
``phi`` (e.g. from a ``--blup`` or ``--no-importance`` estimate) raises an error.

Each row of ``--out`` is ``sample,time,x,y,...`` — the 0-indexed sample (matching the row order in ``--locations``), the time in the past, and one column per spatial dimension of the estimated ancestor location at that time (two columns, ``x`` and ``y``, for the 2D test data). With ``--blup``, each row gets one extra trailing column: the variance of the BLUP estimate.

.. _ancient-samples:

Ancient samples
----------------

As just described, ``process-times`` always reads each sample's age straight back off the
genealogy Relate produced -- there's no separate mode or ages input to ``cli.py`` to opt into.
The walkthrough above happens to use data where ages are all effectively 0; you can run the exact
same five steps on ``data/test_with_ancients.*`` instead -- 
which contains 100 contemporary and 10 ancient diploid individuals.


Predicting a location from relatedness
----------------------------------------

``--forget-locations`` predicts a sample's location from its genetic relatedness to everyone
else alone, under Brownian motion dispersal -- masking the given samples' own known locations
(all held out together) and estimating them purely from the tree::

    python cli.py locate-ancestors \
        --in data/test_chr1_locus1_10000T \
        --locations data/test.locations \
        --sigma data/test.sigma \
        --samples 0 1 \
        --ancestor_times 100 \
        --forget-locations \
        --out data/test_chr1_locus1.forgotten_locs

Two things this is good for: a check of how well the model fits, or a biological question -- e.g.
where a long-distance migrant lineage's ancestors actually came from, based on relatedness rather
than its observed (possibly recently-migrated) location.

Note that asking to locate a sample at ``--ancestor_times`` more recent than a sample's
own age returns ``nan`` (it didn't exist yet), and omitting ``--ancestor_times`` entirely defaults
to each sample's own age, for use with ``--forget-locations``.
