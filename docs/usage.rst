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

See :doc:`api` for the full reference.
