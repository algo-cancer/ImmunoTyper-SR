import argparse, os
from posixpath import splitext
from immunotyper.allele_database import ImgtNovelAlleleDatabase
from immunotyper.common import log, initialize_logger
from immunotyper.bam_filter_classes import IghHg38BamFilter, IghHg37BamFilter
from Bio import SeqIO
from statistics import mean
from immunotyper.read_filter_classes import BwaFlankingFilter
from immunotyper.candidate_builder_classes import BwaMappingsCandidateBuilder
from immunotyper.lpinterface import ShortReadModelTotalErrorDiscardObj
from immunotyper.common import resource_path

parser = argparse.ArgumentParser(description='ImmunoTyper-SR: Ig Genotyping using Short Read WGS')

parser.add_argument('bam_path', type=str, help='Input BAM file')
parser.add_argument('--output_dir', default='', type=str, help='Path to output directory. Outputs txt file of allele calls with prefix matching input BAM file name.')


parser.add_argument('--hg37',
    help='Flag if BAM mapped to GRCh37 not GRCh38',
    action='store_true')

parser.add_argument(
    '--bwa',
    type=str,
    default='bwa',
    help='path to bwa executible if not in $PATH'
)

parser.add_argument(
    '--max_copy',
    type=int,
    default=4,
    help='Maximum number of allele copies to call'
)

parser.add_argument(
    '--landmarks_per_group',
    type=int,
    default=6,
    help='Number of landmarks per group to use (default = 6)'
)

parser.add_argument(
    '--landmark_groups',
    type=int,
    default=6,
    help='Number of landmark groups to use (default = 6)'
)

parser.add_argument(
    '--stdev_coeff',
    type=float,
    default=1.5,
    help='Standard deviation scaling coefficient (default = 1.5)'
)

parser.add_argument(
    '--seq_error_rate',
    type=float,
    default=0.02,
    help='Expected sequence error rate (default = 0.02)'
)

parser.add_argument(
    '--solver_time_limit',
    type=float,
    default=1,
    help='Time limit for ILP solver in hours'
)

parser.add_argument(
    '--debug_log_path',
    default=None,
    help='Path to write log'
)

parser.add_argument(
    '--write_cache_path',
    default=None,
    help='Specific location and name of allele db sam mapping cache'
)

parser.add_argument('--no_coverage_estimation', help='Disables empirical coverage',
    action='store_true')


def main():
    args = parser.parse_args()
    run_immunotyper(args.bam_path, args.hg37, args.output_dir, args.landmark_groups, args.landmarks_per_group, args.max_copy, args.stdev_coeff, args.seq_error_rate, args.write_cache_path, args.solver_time_limit)

def run_immunotyper(bam_path: str,  hg37: bool=False, 
                                    output_dir: str='',
                                    landmark_groups: int=6, 
                                    landmarks_per_group: int=6, 
                                    max_copy: int=4, 
                                    stdev_coeff: float=1.5, 
                                    seq_error_rate: float=0.02,
                                    write_cache_path: str='',
                                    solver_time_limit: int=1):
    allele_db = ImgtNovelAlleleDatabase(**{'db_fasta_path': resource_path('IMGT_IGHV_reference_allele_db-updated+oscar_novel-aligned.fasta'),
                                                                'consensus_path': resource_path('V-QUEST-reference-allele-db+no-period-references+consensus.clustalw.consensus-seq.fasta'),
                                                                'ignored_alleles_path': resource_path('ignored_alleles.txt')})


    # Extract reads from BAM
    if not hg37:
        bam_filter = IghHg38BamFilter(bam_file_path=bam_path)
    else:
        bam_filter = IghHg37BamFilter(bam_file_path=bam_path)
    bam_filter.recruit_reads()
    m, variance, edge_variance = bam_filter.sample_coverage()
    READ_DEPTH = int(round(m/2))
    VARIANCE = variance/2
    EDGE_VARIANCE = [x/2 for x in edge_variance]

    # get read lengths
    reads = SeqIO.parse(bam_filter.output_path, 'fasta')
    lengths = []
    i = 0
    for r in reads:
        lengths.append(len(r))
        i+=1
        if i>100: break
    READ_LENGTH = int(mean(lengths))
    log.info(f"Sampled read length: {READ_LENGTH}")

    allele_db.make_landmarks(landmark_groups*landmarks_per_group, READ_LENGTH, READ_DEPTH, VARIANCE, EDGE_VARIANCE, 50, landmark_groups)

    # Make read mappings to allele database and filter
    flanking_filter = BwaFlankingFilter(reference_path=resource_path('IMGT_IGHV_reference_allele_db-updated+oscar_novel-no_duplicates.fa'),
                                        write_cache_path = write_cache_path if write_cache_path else None,
                                        load_cache_path = write_cache_path if write_cache_path else None)
    filtered_reads, negative = flanking_filter.filter_reads(bam_filter.output_path)


    # Make allele candidates
    candidate_builder = BwaMappingsCandidateBuilder(read_length=READ_LENGTH,
                                                                allele_db=allele_db)
    candidates = candidate_builder.make_candidates(filtered_reads)


    # Build and run model
    model = ShortReadModelTotalErrorDiscardObj(num_landmarks=landmark_groups*landmarks_per_group,
                                            num_landmark_groups=landmark_groups,
                        stdev_coefficient=stdev_coeff, 
                        maxcopy=max_copy,
                            sequencing_error_rate=seq_error_rate)
    model.build(filtered_reads, candidates)
    model.solve(timeLimit=solver_time_limit*3600, threads=12, log_path=None)


    # Write outputs
    output_dir = output_dir if output_dir else os.getcwd()
    output_file = os.path.join(output_dir, os.path.splitext(os.path.basename(bam_path))[0]+'-IGHV_functional_allele_calls.txt')
    log.info(f"Writing allele calls to: {output_file}")
    with open(output_file, 'w') as f:
        for c in model.candidates:
            if allele_db[c.id].is_functional:
                for copy in range(int(c.copy_multipler.getValue())):
                    f.write(c.id+'\n')


if __name__ == '__main__':
	main() 
