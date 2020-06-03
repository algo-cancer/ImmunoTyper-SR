# ImmunoTyper

A Immunoglobulin Heavy Chain Variable Gene genotyping and CNV analysis tool for WGS PacBio Long Reads using ILP Optimization.


## Dependencies:

ImmunoTyper is developed using python 2.7 for compatability with legacy computing clusters.

You can install the following dependencies by running `./install-dependencies.sh`

- [Poa](https://sourceforge.net/projects/poamsa/)
- [SPOA](https://github.com/rvaser/spoa)
- [Dense Subgraph Finder from the IG Repertoire Constructor Toolset](http://yana-safonova.github.io/ig_repertoire_constructor/)
- [Minimap2](https://github.com/lh3/minimap2)

The following python dependencies should be installed using your prefered method (pip, conda etc)
- [Biopython](https://biopython.org)
- [Gurobi](http://www.gurobi.com) via the gurobipy package
- [Blasr](https://github.com/PacificBiosciences/blasr) (easiest to install [via Conda](https://anaconda.org/bioconda/blasr))

You can create a conda environment containing all necessary python devendecies by running 
`conda env create -f environment.yml`
You can then activate the environment using 
`conda activate immunotyper` 


## Running ImmunoTyper:

ImmunoTyper is run from the `src/immunotyper.py` script via the command `python immunotyper.py`.

```
usage: immunotyper.py [-h] [--dsf DSF] [--poa POA] [--spoa SPOA]
                      [--minimap MINIMAP] [--blasr BLASR]
                      [--max_copy MAX_COPY] [--debug_log_path DEBUG_LOG_PATH]
                      [--no_coverage_estimation] [--skip_extraction]
                      reads_path expected_coverage output_dir

ImmunoTyper: IGHV Genotyping using PacBio Long Reads

positional arguments:
  reads_path            Input FASTA file of PacBio reads
  expected_coverage     Sequencing depth
  output_dir            Path to output directory. The directory must exist. A
                        fasta is created for each allele call containing the
                        subreads assigned to that call

optional arguments:
  -h, --help            show this help message and exit
  --dsf DSF             Path to Dense Subgraph Finder tool (usually located in
                        ig_repertoire_constructor/dense_subgraph_finder.py)
  --poa POA             Path to POA tool (usually located at poaV2/poa).
                        Default location is ../tools/poaV2/poa.
  --spoa SPOA           Path to SPOA tool (usually located at
                        spoa/build/bin/spoa). Default location is
                        ../tools/spoa/build/bin/spoa
  --minimap MINIMAP     Path to minimap tool
  --blasr BLASR         Path to blasr tool
  --max_copy MAX_COPY   Maximum number of allele copies to call
  --debug_log_path DEBUG_LOG_PATH
                        Path to write log
  --no_coverage_estimation
                        Disables empirical coverage
  --skip_extraction     Disables extractinbg subreads from input FASTA reads.
                        Used when subreads have already been isolated (i.e. in
                        simulated data)

```
