# spacetrees
Code to estimate dispersal rates and locate genetic ancestors from genome-wide genealogies. Associated with the paper, Osmond & Coop 2024: https://elifesciences.org/articles/72177.

This version of the code has been updated to allow for ancient samples.

Formerly referred to as sparg, but that name is now reserved for inferring spatial histories from full ancestral recombination graphs (https://github.com/osmond-lab/sparg).

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
	- move back to the main working directory, `cd -`
- Install Relate. See https://myersgroup.github.io/relate/index.html for more info and options.
	- On my server I downloaded the source code with `git clone https://github.com/MyersGroup/relate.git` in June 2025. This is roughly version 1.2.2.
	- Move into the Relate directory, `cd relate/build`.
        - On my server I had to load these tools to build relate, `module load cmake/3.22.5 gcc/11.3.0`.
    	- `cmake ..`
    	- `make` 
	- move back to the main working directory, `cd -`
- Run spacetrees via snakemake
	- you should now be able to estimate dispersal and locate genetic ancestors with spacetrees via snakemake! simply write `snakemake all -c1` in the command line (-c1 indicates 1 thread, use more if you have them, but this example should run in less than a minute or two with -c1)
	- on our server we have to first enter `export XDG_CACHE_HOME=$SCRATCH` so that snakemake writes to a write-able location
	- TODO: lots more detail needed about how to customize your options within Snakefile
	- TODO: show how to do this when you have multiple chromosomes
- Plot
	- make virtual environment accessible in Jupyter notebook with `python -m ipykernel install --name $myenv --user` and `venv2jup`
	- TODO: some may need to install Jupyter?
	- open the Jupyter notebook plots.ipynb. I do this through my server's JupyterHub, https://ondemand.scinet.utoronto.ca
	- run the code (command+enter to execute a cell)
- Get in touch!
	- I'd love to hear if you are using this software, have any suggested improvements, or need any help: mm.osmond@utoronto.ca
