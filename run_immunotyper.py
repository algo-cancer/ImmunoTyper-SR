import argparse
from src.immunotyper import igh
from src.immunotyper.common import fasta_from_seq

parser = argparse.ArgumentParser(description='ImmunoTyper: IGHV Genotyping using PacBio Long Reads')
parser.add_argument('reads_path', type=str, help='Input FASTA file of PacBio reads')
parser.add_argument('expected_coverage', type=int, help='Sequencing depth')
parser.add_argument('output_dir', type=str, help='Path to output directory. The directory must exist. A fasta is created for each allele call containing the subreads assigned to that call')


parser.add_argument(
    '--dsf',
    type=str,
    default='dsf',
    help='Path to Dense Subgraph Finder tool (usually located in ig_repertoire_constructor/dense_subgraph_finder.py)'
)
parser.add_argument(
    '--poa',
    type=str,
    default='../tools/poaV2/poa',
    help='Path to POA tool (usually located at poaV2/poa). Default location is ../tools/poaV2/poa.'
)
parser.add_argument(
    '--spoa',
    type=str,
    default='../tools/spoa/build/bin/spoa',
    help='Path to SPOA tool (usually located at spoa/build/bin/spoa). Default location is ../tools/spoa/build/bin/spoa'
)
parser.add_argument(
    '--minimap',
    type=str,
    default='minimap',
    help='Path to minimap tool'
)
parser.add_argument(
    '--blasr',
    type=str,
    default='blasr',
    help='Path to blasr tool'
)

parser.add_argument(
    '--max_copy',
    type=int,
    default=4,
    help='Maximum number of allele copies to call'
)

parser.add_argument(
    '--debug_log_path',
    default=None,
    help='Path to write log'
)

parser.add_argument('--no_coverage_estimation', help='Disables empirical coverage',
    action='store_true')
parser.add_argument('--skip_extraction', help='Disables extractinbg subreads from input FASTA reads. Used when subreads have already been isolated (i.e. in simulated data)',
    action='store_true')


args = parser.parse_args()


def main():

	clusters = igh.run(reads_path=args.reads_path,
						expected_coverage=args.expected_coverage,
						blasr_src=args.blasr,
						minimap=args.minimap,
						dsf_src=args.dsf,
						no_cov_estimation=args.no_coverage_estimation,
						skip_extraction=args.skip_extraction,
						debug_log_path=args.debug_log_path)

	for c in clusters:
		output = args.output_dir if '/' == args.output_dir[-1] else args.output_dir+'/'
		with open('{}{}_number={}.fasta'.format(output, c.call[0], len(c.call)), 'wb') as f:
			f.write(fasta_from_seq('consensus', c.consensus))
			f.write('\n')
			for read in c:
				f.write(fasta_from_seq(read.name, read.seq))
				f.write('\n')
			f.flush()


if __name__ == '__main__':
	main()
