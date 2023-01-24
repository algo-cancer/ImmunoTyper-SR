# ImmunoTyper-SR

An Immunoglobulin Variable Gene genotyping and CNV analysis tool for WGS short reads using ILP Optimization. See our [paper here](https://www.cell.com/cell-systems/fulltext/S2405-4712(22)00352-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471222003520%3Fshowall%3Dtrue).

Note we now support calling IGLV, IGKV and TRAV (see "Running ImmunoTyper-SR" below).


## Installation

If you already have BWA installed and don't want to create a new environment, just download the latest release binary (see right toolbar) and install with pip:

Then run `pip install <binary.whl>`

However BEST way is to setup an clean environment for the installation first (see below).

### Environment and Dependencies

Installing ImmunoTyper-SR using pip will automatically install the following dependencies:

- [biopython](https://biopython.org/)
- [dill](https://pypi.org/project/dill/)
- [gurobipy](https://www.gurobi.com/documentation/9.5/quickstart_mac/cs_grbpy_the_gurobi_python.html)
- [logbook](https://logbook.readthedocs.io/en/stable/)
- [pysam](https://pysam.readthedocs.io/en/latest/api.html)

In addition to the above, you will need 

1.  [BWA mem](http://bio-bwa.sourceforge.net/bwa.shtml) mapper. We recommend using a new conda environment for the installation, which you can also use to install BWA:

```
conda create -n immunotyper-SR -c bioconda python=3.8 bwa samtools
conda activate immunotyper-SR
pip install <binary.whl>
```

2.  [Gurobi](https://www.gurobi.com/) solver configured with a valid license. Licenses are [free for academic purposes](https://www.gurobi.com/downloads/end-user-license-agreement-academic/).

To check that gurobi is correctly configured, run `gurobi_cl` from a shell.

### Installing from source

If the binary fails to install for whatever reason, you can build the tool from source as follows:

```
conda create -n immunotyper-SR -c bioconda python=3.8 bwa samtools
conda activate immunotyper-SR
git clone git@github.com:algo-cancer/ImmunoTyper-SR.git ./ImmunoTyper-SR
cd ImmunoTyper-SR
python -m pip install --upgrade  build
python -m build
pip install dist/<.tar.gz or .whl build>
```



## Running ImmunoTyper-SR:

After installation with pip, simply use the command `immunotyper-SR`. The only required input is a bam file, and the outputs are the following files, made to the current working directory, where `<prefix>` is the input BAM filename without the extension:
-  `<prefix>-<gene_type>_functional_allele_calls.txt`: list of functional alleles called, one per line. Multiple copies are indicated by multiple lines of the same allele.
-  `<prefix>-<gene_type>_allele_calls.txt`: as above but including pseudogenes.
-  `<prefix>-<gene_type>-extracted.fa`: reads extracted from the BAM used for analysis
-  `<prefix>-<gene_type>-immunotyper-debug.log`: log file

IMPORTANT: If your BAM was mapped to GRCh37 use the `--hg37` flag. 

```
$ immunotyper-SR --help
usage: immunotyper-SR [-h] [--gene_type {ighv,iglv,trav,igkv}] [--output_dir OUTPUT_DIR] [--ref REF] [--hg37] [--bwa BWA] [--max_copy MAX_COPY] [--landmarks_per_group LANDMARKS_PER_GROUP] [--landmark_groups LANDMARK_GROUPS] [--stdev_coeff STDEV_COEFF] [--seq_error_rate SEQ_ERROR_RATE] [--solver_time_limit SOLVER_TIME_LIMIT] [--debug_log_path DEBUG_LOG_PATH]
                      [--write_cache_path WRITE_CACHE_PATH] [--threads THREADS] [--no_coverage_estimation]
                      bam_path

ImmunoTyper-SR: Ig Genotyping using Short Read WGS

positional arguments:
  bam_path              Input BAM file

optional arguments:
  -h, --help            show this help message and exit
  --gene_type {ighv,iglv,trav,igkv}
                        Specify which genes to target
  --output_dir OUTPUT_DIR
                        Path to output directory. Outputs txt file of allele calls with prefix matching input BAM file name.
  --ref REF             Path to the reference FASTA to decode CRAM files. Option is not used if bam_path is not a CRAM.
  --hg37                Flag if BAM mapped to GRCh37 not GRCh38
  --bwa BWA             path to bwa executible if not in $PATH
  --max_copy MAX_COPY   Maximum number of allele copies to call
  --landmarks_per_group LANDMARKS_PER_GROUP
                        Number of landmarks per group to use (default = 6)
  --landmark_groups LANDMARK_GROUPS
                        Number of landmark groups to use (default = 6)
  --stdev_coeff STDEV_COEFF
                        Standard deviation scaling coefficient (default = 1.5)
  --seq_error_rate SEQ_ERROR_RATE
                        Expected sequence error rate (default = 0.02)
  --solver_time_limit SOLVER_TIME_LIMIT
                        Time limit for ILP solver in hours
  --debug_log_path DEBUG_LOG_PATH
                        Path to write log
  --write_cache_path WRITE_CACHE_PATH
                        Specific location and name of allele db sam mapping cache
  --threads THREADS     Max number of threads to use
  --no_coverage_estimation
                        Disables empirical coverage
```
