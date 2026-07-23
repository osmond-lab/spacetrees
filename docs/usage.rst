Usage
=====

Running the pipeline
---------------------

``spacetrees`` is driven through `Snakemake <https://snakemake.readthedocs.io/>`_, using the
rules defined in the repository's ``Snakefile``. With the environment from
:doc:`installation` activated, run::

    snakemake all -c1

``-c1`` tells Snakemake to use a single thread; increase this if you have more cores available.
The bundled test data (in ``data/``) should run in well under a minute.

Customizing the run
--------------------

The ``Snakefile`` controls which input trees are used, how dispersal is estimated, and which
ancestor locations are inferred. See the rules and parameters at the top of ``Snakefile`` in
the repository root to adjust these for your own data.

.. note::
   Fuller documentation of the available ``Snakefile`` options, and an example with multiple
   chromosomes, is planned but not yet written — see the repository's README for the latest
   status.

Running core functions from the command line
---------------------------------------------

``cli.py`` exposes :func:`spacetrees.estimate_dispersal` and :func:`spacetrees.locate_ancestors`
directly from the command line, as an alternative to driving them through Snakemake — useful for
scripting or for rerunning inference on already-preprocessed data. It expects the per-locus
``.stss``, ``.stss_logdet``, ``_stss_inv.npy``, ``.btss``, and ``.lpcs`` files produced by the
``Snakefile``'s ``process_times`` rule, plus a sample locations file.

Estimate a dispersal rate from one or more loci::

    python cli.py estimate-dispersal \
        --stss-logdet locus1.stss_logdet locus2.stss_logdet \
        --stss-inv locus1_stss_inv.npy locus2_stss_inv.npy \
        --btss locus1.btss locus2.btss \
        --lpcs locus1.lpcs locus2.lpcs \
        --locations test.locations \
        --out test.sigma

Locate ancestors at a locus, given that estimated rate::

    python cli.py locate-ancestors \
        --stss locus1.stss \
        --stss-inv locus1_stss_inv.npy \
        --btss locus1.btss \
        --lpcs locus1.lpcs \
        --locations test.locations \
        --sigma test.sigma \
        --samples 0 1 \
        --times 10 100 1000 \
        --out locus1.locs

Add ``--blup`` (and optionally ``--blup-var``) to ``locate-ancestors`` for the faster best linear
unbiased predictor instead of the full likelihood surface, matching the Snakefile's
``locate_ancestors_blup`` rule. Run either subcommand with ``--help`` for the full list of options.

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

``cli.py`` (see above) wraps the two ``spacetrees`` functions for direct command-line use.
See :doc:`api` for the full reference.
