Installation
============

These steps set up ``spacetrees`` and its two external dependencies (tsconvert and Relate)
from the command line.

1. Clone the repository and move into it::

    git clone https://github.com/osmond-lab/spacetrees.git
    cd spacetrees

2. Install Python 3.11.5 (or a similar 3.11.x version). On the author's server this is
   done with::

    module load NiaEnv/2022a python/3.11.5

3. Create and activate a virtual environment with that Python version::

    python -m venv venv
    source venv/bin/activate

4. Install the Python dependencies::

    pip install -r requirements.txt

5. Install `tsconvert <https://github.com/tskit-dev/tsconvert>`_, which isn't available on PyPI::

    git clone https://github.com/tskit-dev/tsconvert.git
    cd tsconvert
    pip install .
    cd -

6. Install `Relate <https://myersgroup.github.io/relate/index.html>`_ (v1.2.1 or similar)::

    git clone https://github.com/MyersGroup/relate.git
    cd relate/build
    # module load cmake/3.22.5 gcc/11.3.0   # if needed to build Relate
    cmake ..
    make
    cd -

7. (Optional, for plotting) Make the virtual environment available as a Jupyter kernel::

    python -m ipykernel install --name $myenv --user

Once these are installed you're ready to run ``spacetrees`` — see :doc:`usage`.
