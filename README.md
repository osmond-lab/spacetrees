# spacetrees
Code to estimate dispersal rates and locate genetic ancestors from genome-wide genealogies

# set up

Here is how to get set up and run spacetrees, from the command-line:

- Clone this directory, `git clone https://github.com/osmond-lab/spacetrees.git`.
- Move into this directory, `cd spacetrees`.
- Install Python v3.11.5 (https://www.python.org/downloads/release/python-3115/). On my server we can do this with `module load NiaEnv/2022a python/3.11.5`. May also work with similar versions. 
- Create virtual environment, `python -m venv venv`. Make sure you are using the correct version of Python to do this.
- Activate virtual environment, `source venv/bin/activate`.
- Install Python packages, `pip install -r requirements.txt`.
- Install tsconvert, which isn't available via pip.
	- `git clone https://github.com/tskit-dev/tsconvert.git`. This was v0.1.dev57+g057435c for me, June 7, 2024.
	- `cd tsconvert`
	- `pip install .`
- Install Relate. I used v1.2.1. See https://myersgroup.github.io/relate/index.html for more info and options.
	- On my server I downloaded the source code with `git clone https://github.com/MyersGroup/relate.git` on June 7, 2024. This is roughly version 1.2.1.
	- Move into the Relate directory, `cd relate/build`.
        - On my server I had to load these tools to build relate, `module load cmake/3.22.5 gcc/11.3.0`.
    	- `cmake ..`
    	- `make` 

