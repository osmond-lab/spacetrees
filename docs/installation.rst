Installation
============

These steps set up ``spacetrees`` from the command line.

1. Clone the repository and move into it::

    git clone https://github.com/osmond-lab/spacetrees.git
    cd spacetrees

2. Create and activate a virtual environment with python3::

    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt

This also installs `tsconvert <https://github.com/tskit-dev/tsconvert>`_ directly from GitHub (it isn't
on PyPI, so ``requirements.txt`` points pip at a pinned commit instead).

Once this is installed you're ready to run ``spacetrees`` — see :doc:`usage`.
