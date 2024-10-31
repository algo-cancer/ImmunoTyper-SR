import argparse, os
from posixpath import splitext
from .common import log, initialize_logger, get_database_config, get_allele_db_mapping_path
from .allele_database import ImgtNovelAlleleDatabase
from .bam_filter_classes import BamFilterImplemented, IghHg38BamFilter
from Bio import SeqIO
from statistics import mean
from .read_filter_classes import BowtieFlankingFilter
from .candidate_builder_classes import BwaMappingsCandidateBuilder
from .solvers import GurobiSolver
from .models import ShortReadModelTotalErrorDiscardObj
from .common import resource_path
from .post_processing import PostProcessorModel


implemented_solvers = {'gurobi': GurobiSolver}

implemented_models = {'ilp': ShortReadModelTotalErrorDiscardObj}


parser = argparse.ArgumentParser(description='ImmunoTyper-SR: Ig Genotyping using Short Read WGS')

parser.add_argument('bam_path', type=str, help='Input BAM file')
parser.add_argument('--gene_type', choices=['ighv', 'iglv', 'trav', 'igkv', 'trbv', 'trdv', 'trgv'], default='ighv', help='Specify which genes to target')
parser.add_argument('--output_dir', default='', type=str, help='Path to output directory. Outputs txt file of allele calls with prefix matching input BAM file name.')


parser.add_argument(
    '--ref',
    type=str,
    help='Path to the reference FASTA to decode CRAM files. Option is not used if bam_path is not a CRAM.',
)

parser.add_argument('--hg37',
    help='Flag if BAM mapped to GRCh37 not GRCh38',
    action='store_true')

parser.add_argument('--solver',
    help='Choose ilp solver',
    choices=['gurobi'], 
    default='gurobi'
)

# parser.add_argument('--model',
#     help='Choose model type',
#     choices=['ilp'], 
#     default='ilp'
# )

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
    type=str,
    default=None,
    help='Specific location and name of allele db sam mapping cache'
)

parser.add_argument(
    '--threads',
    type=int,
    default=6,
    help='Max number of threads to use'
)


parser.add_argument('--no_coverage_estimation', help='Disables empirical coverage',
    action='store_true')

parser.add_argument(
    '--save_extracted_reads',
    action='store_true',
    help='Save the extracted reads FASTA file instead of deleting it after use'
)

def main():
    args = parser.parse_args()    
    run_immunotyper(args.bam_path, args.ref, args.gene_type, args.hg37, args.solver, args.output_dir, args.landmark_groups, args.landmarks_per_group, args.max_copy, args.stdev_coeff, args.seq_error_rate, args.write_cache_path, args.solver_time_limit, args.threads, args.save_extracted_reads)


def run_immunotyper(bam_path: str,  ref: str='',
                                    gene_type: str='ighv',
                                    hg37: bool=False, 
                                    solver: str='gurobi',
                                    output_dir: str='',
                                    landmark_groups: int=6, 
                                    landmarks_per_group: int=6, 
                                    max_copy: int=4, 
                                    stdev_coeff: float=1.5, 
                                    seq_error_rate: float=0.02,
                                    write_cache_path: str='',
                                    solver_time_limit: int=1,
                                    threads: int=6,
                                    save_extracted_reads: bool=False):
    """Driver method to run immunotyper and output calls

    Args:
        bam_path (str): Path to Input BAM file
        gene_type (str, optional): Defaults to 'ighv'.
        hg37 (bool, optional): Flag if BAM mapped to GRCh37. Defaults to False.
        model_type (str, optional): Model to use. Defaults to 'ilp-gurobi'.
        output_dir (str, optional): _description_. Defaults to ''.
        landmark_groups (int, optional): _description_. Defaults to 6.
        landmarks_per_group (int, optional): _description_. Defaults to 6.
        max_copy (int, optional): _description_. Defaults to 4.
        stdev_coeff (float, optional): _description_. Defaults to 1.5.
        seq_error_rate (float, optional): _description_. Defaults to 0.02.
        write_cache_path (str, optional): _description_. Defaults to ''.
        solver_time_limit (int, optional): _description_. Defaults to 1.
        threads (int, optional): _description_. Defaults to 6.
    """
    output_prefix = os.path.splitext(os.path.basename(bam_path))[0]
    initialize_logger(os.path.join(output_dir, f'{output_prefix}-{gene_type}-immunotyper-debug'))

    allele_db = ImgtNovelAlleleDatabase(**get_database_config(gene_type))

    # Extract reads from BAM
    bam_filter = BamFilterImplemented(bam_path, gene_type, not hg37, reference_fasta_path=ref, output_path=output_dir)
    reads = bam_filter.recruit_reads()
    m, variance, edge_variance = bam_filter.sample_coverage(large_depth_sample=True)
    READ_DEPTH = int(round(m/2))
    VARIANCE = variance/2
    EDGE_VARIANCE = [x/2 for x in edge_variance]

    # sample read lengths
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
    flanking_filter = BowtieFlankingFilter(reference_path=get_allele_db_mapping_path(gene_type),
                                    write_cache_path = write_cache_path if write_cache_path else None,
                                    load_cache_path = write_cache_path if write_cache_path else None)
                                   
    positive, negative = flanking_filter.filter_reads(bam_filter.output_path, mapping_params="-a --end-to-end --very-sensitive -f  --n-ceil C,100,0 --np 0 --ignore-quals --mp 2,2 --score-min C,-50,0 -L 10")
    for r in positive: r.allele_db = allele_db
    
    if not save_extracted_reads:
        bam_filter.delete_extracted_reads()
    
    # Make allele candidates
    ### Instantiate model to get candidate class to use
    candidate_builder = BwaMappingsCandidateBuilder(read_length=READ_LENGTH,
                                                                allele_db=allele_db)
    candidates = candidate_builder.make_candidates(positive)


    # Build and run model
    model = implemented_models['ilp'](implemented_solvers[solver], num_landmarks=landmark_groups*landmarks_per_group,
                                            num_landmark_groups=landmark_groups,
                        stdev_coefficient=stdev_coeff, 
                        maxcopy=max_copy,
                            sequencing_error_rate=seq_error_rate)

    model.build(positive, candidates)
    model.solve(time_limit=solver_time_limit*3600, threads=threads, log_path=os.path.join(output_dir, f'{output_prefix}-{gene_type}-gurobi.log'))


    # Write outputs
    output_dir = output_dir if output_dir else os.getcwd()
    
    # Write functional allele calls
    output_file = os.path.join(output_dir, os.path.splitext(os.path.basename(bam_path))[0]+f'-{gene_type.upper()}_functional_allele_calls.txt')
    log.info(f"Writing functional allele calls to: {output_file}")
    model.write_allele_calls(output_file, functional_only=True)
    
    # Write all allele calls
    output_file = os.path.join(output_dir, os.path.splitext(os.path.basename(bam_path))[0]+f'-{gene_type.upper()}_allele_calls.txt')
    log.info(f"Writing all allele calls to: {output_file}")
    model.write_allele_calls(output_file, functional_only=False)
            
    # Call SNVs
    post_processor = PostProcessorModel(model, allele_db, sequencing_depth=READ_DEPTH)

    # Ensure the output VCF directory exists
    output_vcf_dir = os.path.join(output_dir, os.path.splitext(os.path.basename(bam_path))[0]+'-novel_variant_vcfs')
    os.makedirs(output_vcf_dir, exist_ok=True)

    output_path = os.path.join(output_dir, os.path.splitext(os.path.basename(bam_path))[0]+'-novel-variants.tsv')
    log.info(f'Writing novel variants to: {output_path}')
    
    post_processor.write_all_variants(
        output_path=output_path,
        output_vcf_dir=output_vcf_dir
    )

    # Write read assignments
    output_dir = os.path.join(output_dir, os.path.splitext(os.path.basename(bam_path))[0]+'-read_assignment')
    log.info(f'Writing read assignments to: {output_dir}')
    post_processor.write_read_assignments(
        output_dir=output_dir
    )

if __name__ == '__main__':
	main() 

