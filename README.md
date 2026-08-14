# spacetrees
Code to estimate dispersal rates and locate genetic ancestors from genome-wide genealogies. Associated with the paper, Osmond & Coop 2024: https://elifesciences.org/articles/72177.

Formerly referred to as sparg, but that name is now reserved for inferring spatial histories from full ancestral recombination graphs (https://github.com/osmond-lab/sparg).

# warning
If you need a time cutoff use the manuscript code (https://github.com/mmosmond/spacetrees-ms), we need to fix a bug in this streamlined version.

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
- Install Relate. Go to https://myersgroup.github.io/relate/ and download a precompiled binary (pick the static build unless you have a specific reason to prefer dynamic — it doesn't depend on your system's library versions matching what Relate was built with).
	- Extract the downloaded archive into the project directory and rename the resulting folder to `relate` (i.e. so `relate/scripts/...` and `relate/bin/...` exist at the project root).
	- If you'd rather build from source, or need a version not offered as a precompiled binary, see https://github.com/MyersGroup/relate for instructions.
- Run spacetrees from the command line
	- `cli.py` covers the whole pipeline from a Relate `.mut` file to located ancestors, as five subcommands: `python cli.py loci-positions --mut test_chr1.mut --out test_chr1.loci`, then (once you've sampled trees at a locus with `SampleBranchLengths.sh`) `python cli.py extract-times --newick locus1.newick --out locus1` (writes `locus1.stss` and `locus1.ctss`), then `python cli.py process-times --times locus1 --coal test.coal --T 10000 --out locus1` (reads `locus1.stss` and `locus1.ctss`; `--out` is a base prefix, so this writes `locus1_10000T.stss`, `locus1_10000T.stss-logdet`, `locus1_10000T_stss-inv.npy`, `locus1_10000T.btss`, and `locus1_10000T.lpcs` — omit `--T` to skip the `.stss` output, which would just duplicate `extract-times`', and write the rest straight to `locus1.stss-logdet`/`locus1_stss-inv.npy`/`locus1.btss`/`locus1.lpcs`, no suffix)
	- once you have those preprocessed per-locus files, `python cli.py estimate-dispersal --in locus1_10000T locus2_10000T --locations test.locations --out test.sigma` (reads `locus1_10000T.stss-logdet`/`locus1_10000T_stss-inv.npy`/`locus1_10000T.btss`/`locus1_10000T.lpcs` and the same for `locus2_10000T`) estimates a dispersal rate, and `python cli.py locate-ancestors --in locus1_10000T --locations test.locations --sigma test.sigma --ancestor_times 10 100 1000 --out locus1.locs` locates ancestors. Run `python cli.py <subcommand> --help` for the full set of options for any of the five subcommands (`loci-positions`, `extract-times`, `process-times`, `estimate-dispersal`, `locate-ancestors`).
- Plot
	- make virtual environment accessible in Jupyter notebook with `python -m ipykernel install --name $myenv --user` and `venv2jup`
	- TODO: some may need to install Jupyter?
	- open the Jupyter notebook plots.ipynb. I do this through my server's JupyterHub, https://jupyter.scinet.utoronto.ca/
	- run the code (command+enter to execute a cell)
- Get in touch!
	- I'd love to hear if you are using this software, have any suggested improvements, or need any help: mm.osmond@utoronto.ca
