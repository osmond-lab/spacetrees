Usage
=====

Running the pipeline
---------------------

``spacetrees`` is driven through ``cli.py``, a command-line interface with five subcommands
covering the whole pipeline from a Relate ``.mut`` file to located ancestors: ``loci-positions``,
``extract-times``, ``process-times``, ``estimate-dispersal``, and ``locate-ancestors``. With the
environment from :doc:`installation` activated, run any of them as::

    python cli.py <subcommand> --help

to see its full list of options. You can jump in wherever you already have intermediate files —
each subcommand only needs the files described below (the ``sample_trees`` step itself, i.e.
Relate's ``SampleBranchLengths.sh``, is not wrapped here — run it directly to produce the newick
file each locus's trees are sampled into).

Get the position of the first and last mutation at each locus in a ``.mut`` file::

    python cli.py loci-positions \
        --mut test_chr1.mut \
        --out test_chr1.loci

Extract shared and coalescence times from the sampled trees at a locus::

    python cli.py extract-times \
        --newick locus1.newick \
        --out locus1

This writes ``locus1.stss`` and ``locus1.ctss``.

Process those times at a given time cutoff ``T`` (chop, center, invert the shared times, and
derive branching times and coalescent-time log probabilities; omit ``--T`` for no cutoff).
``--out`` is a base prefix — the cutoff is appended automatically, so the files this writes
depend on ``--T``::

    python cli.py process-times \
        --times locus1 \
        --coal test.coal \
        --T 10000 \
        --out locus1

This reads ``locus1.stss`` and ``locus1.ctss`` (from ``extract-times`` above) and writes
``locus1_10000T.stss``, ``locus1_10000T.stss-logdet``, ``locus1_10000T_stss-inv.npy``,
``locus1_10000T.btss``, and ``locus1_10000T.lpcs`` — everything under one self-contained
``locus1_10000T`` prefix, ready for :func:`spacetrees.estimate_dispersal` and
:func:`spacetrees.locate_ancestors` below (plus a sample locations file). Different ``--T``
values with the same ``--out`` don't collide — each gets its own suffixed prefix.

Omitting ``--T`` skips the ``.stss`` output — chopping is a no-op with no cutoff, so it would only
duplicate ``extract-times``' file — and writes the rest directly to ``--out`` with no suffix. Use
the *same* prefix as ``--times`` so the unchopped ``locus1.stss`` sits right alongside them::

    python cli.py process-times \
        --times locus1 \
        --coal test.coal \
        --out locus1

This writes ``locus1.stss-logdet``, ``locus1_stss-inv.npy``, ``locus1.btss``, and ``locus1.lpcs``
— combined with the existing ``locus1.stss``, ``--in locus1`` below has everything it needs.

Estimate a dispersal rate from one or more loci::

    python cli.py estimate-dispersal \
        --in locus1_10000T locus2_10000T \
        --locations test.locations \
        --out test.sigma

Each ``--in`` entry is a full prefix from ``process-times`` (including its ``T`` suffix): this
reads ``locus1_10000T.stss-logdet``, ``locus1_10000T_stss-inv.npy``, ``locus1_10000T.btss``, and
``locus1_10000T.lpcs``, and the same four files for ``locus2_10000T``.

Locate ancestors at a locus, given that estimated rate::

    python cli.py locate-ancestors \
        --in locus1_10000T \
        --locations test.locations \
        --sigma test.sigma \
        --samples 0 1 \
        --ancestor_times 10 100 1000 \
        --out locus1.locs

``--in`` reads all four files written by ``process-times`` for that prefix: ``locus1_10000T.stss``
(already chopped at ``T=10000``), ``locus1_10000T_stss-inv.npy``, ``locus1_10000T.btss``, and
``locus1_10000T.lpcs``.

Add ``--blup`` (and optionally ``--blup-var``) to ``locate-ancestors`` for the faster best linear
unbiased predictor instead of the full likelihood surface. Run any subcommand with ``--help`` for
its full list of options.

Plotting results
-----------------

Results are explored in the bundled ``plots.ipynb`` notebook:

1. Make the project's virtual environment available to Jupyter (see :doc:`installation`).
2. Open ``plots.ipynb`` (e.g. through JupyterHub) and run the cells to visualize dispersal
   estimates and located ancestors.

Core functions
--------------

The pipeline calls into two Python modules for the underlying computation:

* :mod:`spacetrees` — :func:`spacetrees.estimate_dispersal` fits the maximum likelihood
  dispersal (and optionally branching) rate from sample locations and shared coalescence
  times; :func:`spacetrees.locate_ancestors` uses a fitted dispersal rate to locate genetic
  ancestors at given times.
* :mod:`utils` — helpers for extracting shared coalescence times from trees and computing
  coalescent time densities.

``cli.py`` (see above) wraps ``utils``' tree-processing helpers and the two ``spacetrees``
functions for direct command-line use. See :doc:`api` for the full reference.
