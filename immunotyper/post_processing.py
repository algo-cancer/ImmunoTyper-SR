import os, pysam
from tempfile import TemporaryFile
from typing import Tuple, List
from venv import create
from .common import create_temp_file as ctf
from .common import fasta_from_seq, suppress_stdout, Read, log, run_command
from .mapper_wrappers import BowtieWrapper
from statistics import mean
from Bio import SeqIO
from glob import glob

def create_temp_file(*args, **kwargs):
    kwargs['delete'] = False
    return ctf(*args, **kwargs)


class PostProcessor(object):

    mapper = BowtieWrapper(params='--very-sensitive-local')

    # gene_mapped_reads_bam_paths = {}        # path names of bam files containing mappings of all reads assigned to gene
    # reference_target_fa_paths = {}          # path names of indexed mapping target allele sequences
    # vcf_variants = {}
    gene_type = 'IGHV'


    def __init__(self, sequencing_depth: int, read_assignment, input_reads, calls, allele_db, minimum_coding_bases: int = 50) -> None:
        '''
        Args
            sequencing_depth            Single copy / haploid sequencing depth'''
        self.gene_mapped_reads_bam_paths = {}
        self.reference_target_fa_paths = {}
        self.vcf_variants = {}
        self.sequencing_depth = sequencing_depth
        self.minimum_coding_bases = minimum_coding_bases
        self.read_assignment = read_assignment
        self.input_reads = input_reads
        self.allele_db = allele_db
        self.calls = calls
        self.set_calls(self.calls)


    def set_calls(self, raw_calls):
        self.calls_no_cnv = list(set(raw_calls))
        self.gene_calls = list(set([x.split('*')[0] for x in raw_calls]))
        self.functional_calls = [c for c in raw_calls if self.allele_db[c].is_functional]
        self.functional_calls_no_cnv = list(set(self.functional_calls))
        self.functional_gene_calls = list([x.split('*')[0] for x in self.functional_calls])

    def check_call_has_snv_support(self, gene: str, problem_allele: str, mapping_target: str, required_snps: list, override_gene_alleles: List[str] = []) -> bool:
        ''' Checks if set of distinguishing variants for problem_allele are found in the set of variants called from 
        all reads assigned to alleles of problem_allele's parent gene
        Args
            gene (str)                          Gene with which to collect assigned reads
            problem_allele (str)                Allele to be replace if required_snps not present in vcf
            mapping_target (str)                Target to map reads against for variant calling. Usually wildtype
            required_snps (tuple(pos, var_val)) Variants required to be present to NOT change problem_allele into replacement_allele
            override_gene_alleles (list)        Alleles to collect reads from. If == None, collects reads from call alleles of gene
            '''
        # check if problem allele in called allele set
        alleles = override_gene_alleles if override_gene_alleles else [a for a in self.calls if gene == self.allele_db[a].gene]
        if problem_allele not in alleles:
            raise ValueError(f"Checking for snv support for {problem_allele}, not present in called alleles: {str(alleles)}")

        try:
            variants = self.vcf_variants[mapping_target][frozenset(alleles)]
        except KeyError:
            self.call_snvs(gene, mapping_target, override_gene_alleles)
            variants = self.vcf_variants[mapping_target][frozenset(alleles)]            

        # Check required_snps are present
        if not all([(pos in variants and var_val in variants[pos]) for pos, var_val in required_snps]): # not all required_snps are present
            return False
        else:
            return True
        
    def call_snvs(self, gene: str, mapping_target: str, override_gene_alleles: List[str] = []):
        '''Calls variants against mapping_target, using reads from called alleles of gene, or alleles in override_gene_alleles
        Sets self.vcf_variants'''
        alleles = override_gene_alleles if override_gene_alleles else [a for a in self.calls if gene == self.allele_db[a].gene]

        # Check if mapping already performed, otherwise call self.map_gene_assignments_to_wildtype
        try:
            mapping_target_path = self.reference_target_fa_paths[mapping_target]
            bam_path = self.gene_mapped_reads_bam_paths[mapping_target][frozenset(alleles)]
        except KeyError:
            self.map_gene_assignments_to_wildtype(gene, mapping_target, override_gene_alleles)
            mapping_target_path = self.reference_target_fa_paths[mapping_target]
            bam_path = self.gene_mapped_reads_bam_paths[mapping_target][frozenset(alleles)]

        # make vcf
        bcf = create_temp_file(suffix='.bcf')
        bcf_path = ''.join([x for x in bcf.name])
        command = f"bcftools mpileup -f {mapping_target_path} {bam_path} | bcftools call -mv -Ob -o {bcf_path}"
        run_command(command, check=True)

        try:
            self.vcf_variants[mapping_target][frozenset(alleles)] = dict([(x.pos, (x.ref, x.alts)) for x in pysam.VariantFile(str(bcf_path))])
        except KeyError:
            self.vcf_variants[mapping_target] = {frozenset(alleles): dict([(x.pos, (x.ref, x.alts)) for x in pysam.VariantFile(str(bcf_path))])}
        except ValueError as e:
            log.error(f"Error running {command}:")
            raise e

    
    def call_all_snvs(self):
        for gene in self.functional_gene_calls:
            if gene not in self.allele_db.gene_in_cluster:
                mapping_target = self.allele_db.get_wildtype_target(gene)
                # print(f"Calling variants for {gene} against {mapping_target}")
                try:
                    self.call_snvs(gene, mapping_target)
                except Exception as e:
                    log.warn(f"Exception calling SNVs for {gene} against {mapping_target}: {str(e)}")
                    raise e
        for cluster in self.allele_db.gene_clusters:
            called = [g for g in cluster if g in self.gene_calls]
            if called:
                alleles = [a for a in self.functional_calls if self.allele_db[a].gene in called]
                mapping_target = self.allele_db.get_wildtype_target(cluster[0])
                try:
                    self.call_snvs(called[0], mapping_target, override_gene_alleles=alleles)
                except Exception as e:
                    log.warn(f"Exception calling SNVs for {str(called)} group against {mapping_target}: {str(e)}")



    def map_gene_assignments_to_wildtype(self, gene: str, mapping_target: str, override_gene_alleles: List[str] = []) -> None:
        '''
        Args
            gene (str)                          Gene with which to collect assigned reads
            mapping_target (str)                Target to map reads against for variant calling. Usually wildtype
            override_gene_alleles (list)        Alleles to collect reads from. If == None, collects reads from call alleles of gene
            '''
        
        # collect reads
        if override_gene_alleles:
            gene_assigned_reads = [r for a in override_gene_alleles for r in self.read_assignment[a]]
            alleles = override_gene_alleles
        else: # get reads from all called alleles of gene
            alleles = [a for a in self.calls if gene == self.allele_db[a].gene]
            gene_assigned_reads = [r for a in alleles for r in self.read_assignment[a]]
        # print(f"Alleles = {str(alleles)}")
        
        _, reference_fa, bam_path = self.map(gene_assigned_reads, self.allele_db[mapping_target])

        try:
            self.gene_mapped_reads_bam_paths[mapping_target][frozenset(alleles)] = bam_path
        except KeyError:
            self.gene_mapped_reads_bam_paths[mapping_target] = {frozenset(alleles): bam_path}
        self.reference_target_fa_paths[mapping_target] = reference_fa

    def map(self, query: List[Read], target: Read) -> Tuple[str, str, str]:
        query = create_temp_file(write_data=fasta_from_seq(*zip(*[(x.id, x.seq) for x in query])), suffix='.fa')
        target = create_temp_file(write_data=fasta_from_seq(target.id, target.seq))
        sam = create_temp_file(suffix='.sam')
        sam_path = sam.name
        
        # map reads, call variants
        with suppress_stdout():
            self.mapper.index_reference(target.name)
            self.mapper.map(query.name, target.name, output_path=sam_path)
        run_command(f"samtools view -b {sam_path} | samtools sort - > {sam_path.replace('sam', 'bam')}")
        sam.close()
        return query.name, target.name, sam_path.replace('sam', 'bam')


    def calculate_called_vs_expected_depth(self, gene, override_gene_alleles: List[str] = []):
        '''
        Args
            override_gene_alleles (list)        Alleles to collect reads from. If == None, collects reads from call alleles of gene
        '''

        called_alleles = override_gene_alleles if override_gene_alleles else [a for a in self.calls if gene == self.allele_db[a].gene]

        if gene+'*01' in self.allele_db:
            mapping_target=gene+'*01'
        else:
            mapping_target = sorted(called_alleles)[0]
        
        try:
            bam_path = self.gene_mapped_reads_bam_paths[mapping_target][frozenset(called_alleles)]
        except KeyError:
            self.map_gene_assignments_to_wildtype(gene, mapping_target, override_gene_alleles)
            bam_path = self.gene_mapped_reads_bam_paths[mapping_target][frozenset(called_alleles)]
        
        # get depth array
        depth = [int(x.split('\t')[2]) for x in pysam.depth(bam_path).strip().split('\n')]

        return len(called_alleles), round(mean(depth[self.minimum_coding_bases:-self.minimum_coding_bases])/self.sequence_depth)


    # def remove_copies(self, gene: str, desired_copy_number: int, alleles_to_remove: list = [], override_gene_alleles: List[str] = []):
    #     '''Removes allele calls from gene to result in desired_copy_number. 
    #     Starts with alleles_to_remove, otherwise removes copies with the highest objectove function error in self.model
    #     Agrs
    #         override_gene_alleles               List of alleles to consider as a "gene"'''
        

    def get_all_variants(self):
        '''Formats all variants from self.vcf_variants into list of strings of the form "gene_id:positions:alt_nucleotide
        Ensures duplicates dont happen by naive method: ensure no allele is includes in a vcf call set more than once, otherwise 
        raise ValueError

        Returns list of variant strings'''
        if not self.vcf_variants:
            raise AttributeError(f"SNVs not called!")
        variant_labels = []
        alleles_processed = set() # for ensuring no allele contributes variants more than once
        for mapping_target, targets in self.vcf_variants.items():
            gene = self.allele_db[mapping_target].gene
            for alleles, var in targets.items():
                alleles = set([x for x in alleles])
                already_processed = alleles.intersection(alleles_processed)
                if already_processed:
                    raise ValueError(f"Mapping target {mapping_target} has >= 1 allele present in some other variant call")
                already_processed = already_processed.union(alleles)
                for pos, (ref, alts) in var.items():
                    for a in alts:
                        variant_labels.append(".".join([gene, str(pos), ref+a]))
        return variant_labels
    
    def write_all_variants(self, output_path):
        '''Writes all SNV variants provided by self.get_all_variants to output_path. Calls self.call_all_snvs if needed'''
        try:
            var = self.get_all_variants()
        except AttributeError:
            self.call_all_snvs()
            var = self.get_all_variants()
        
        with open(output_path, 'w') as f:
            for v in var:
                f.write(v+'\n')
            

class PostProcessorModel(PostProcessor):
    '''Child of PostProcessor using lpinterface.ShortReadModelTotalErrorDiscardObj as input to get required data'''

    def __init__(self, model, allele_db, sequencing_depth: int=None, *args, **kwargs) -> None:
        '''
        Args
            sequencing_depth            Single copy sampled sequencing depth'''
        super_attr = self.get_attributes_from_model(model)
        kwargs.update(super_attr)
        kwargs['allele_db'] = allele_db
        kwargs['sequencing_depth'] = sequencing_depth
        super().__init__(*args, **kwargs)


    def get_attributes_from_model(self, model):
        # set attr not used in parent class
        self.candidates = dict([(c.id, c) for c in model.candidates])

        read_assignment = {}
        for r in model.reads:
            assigned = [allele_id for allele_id, d_var in r.assignment.items() if d_var.X > 0.5]
            if len(assigned) > 1:
                raise ValueError(f'Read {r.id} assigned to > 1 allele: {str(assigned)}')
            r.assigned = assigned[0] if assigned else None
            if r.assigned:
                try:
                    read_assignment[r.assigned].append(r)
                except KeyError:
                    read_assignment[r.assigned] = [r]
        reads = dict([(r.id, r) for r in model.reads])
        
        calls = []
        for c in self.candidates.values():
            calls.extend([c.id]*int(c.copy_multipler.getValue()))
        
        return {'read_assignment': read_assignment, 'input_reads': reads, 'calls': calls}


class PostProcessorOutput(PostProcessor):
    '''Child of PostProcessor using saved read assignment fastas as input'''

    allele_id_parser = staticmethod(lambda fa_path: fa_path.split('/')[-1].replace('.fa', '').replace('-OR', '/OR'))

    def __init__(self, assigned_read_dir: str, calls_path: str, allele_db, *args, **kwargs):
        self.allele_db = allele_db

        input_reads, read_assignment = self.get_reads(assigned_read_dir)

        calls = self.get_calls(calls_path)

        kwargs.update({'read_assignment': read_assignment,
                        'allele_db': allele_db,
                        'input_reads': input_reads,
                        'calls': calls
                        })

        super().__init__(*args, **kwargs)

    def get_reads(self, assigned_read_dir: str):
        '''Iterates through genotype read assignment fastas provided in assigned_read_dir
        Creates input_reads, read_assignment dictionarys for super.__init__'''
        read_assignment = {}
        reads = {}
        for allele_reads_path in glob(os.path.join(assigned_read_dir, '*.fa')):
            if 'discarded-reads' in allele_reads_path:
                continue
            try: # check assigned allele in fasta is valid
                allele_id = self.allele_id_parser(allele_reads_path)
                if allele_id not in self.allele_db: raise IndexError
            except IndexError as e: # Assigned allele not valid, log and skip
                log.warn(f"Could not load assigned reads from {allele_reads_path}")
                continue
            
            for r in SeqIO.parse(allele_reads_path, 'fasta'): # add read sequences
                try:
                    read_assignment[allele_id].append(r)
                except KeyError:
                    read_assignment[allele_id] = [r]
                reads[r.id] = r
        return reads, read_assignment
        

    def get_calls(self, calls_path):
        calls = []
        with open(calls_path, 'r') as f:
            for line in set([x.strip() for x in f.readlines()]):
                if line[:4] == 'Gene': continue # header
                gene, copies, alleles, _ = line.split('\t')
                copies = int(copies)
                alleles = alleles.split(', ')
                if len(alleles) == copies:
                    calls.extend([self.gene_type+gene+a for a in alleles])
                elif len(alleles) == 1:
                    calls.extend([self.gene_type+gene+alleles[0]]*copies)
                else:
                    raise ValueError(f"Gene call for {gene} in {calls_path} has {copies} copies but lists {', '.join(alleles)} alleles")
        for a in calls:
            if a not in self.allele_db:
                raise ValueError(f"Allele call {a} from {calls_path} not in database")
        return calls
