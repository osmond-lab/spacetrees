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

3. Install `tsconvert <https://github.com/tskit-dev/tsconvert>`_, which isn't available on PyPI::

    git clone https://github.com/tskit-dev/tsconvert.git
    cd tsconvert
    pip install .
    cd -

Once these are installed you're ready to run ``spacetrees`` — see :doc:`usage`.
