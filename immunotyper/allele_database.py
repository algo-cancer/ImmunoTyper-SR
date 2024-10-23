import os, math, json
from Bio import SeqIO
from Bio.Seq import Seq
from abc import ABC, abstractmethod
from .common import log
from itertools import combinations
from collections import OrderedDict, Iterator, defaultdict
from math import sqrt
import dill as pickle
#
# Database classes
#

class AlleleDatabase(ABC, Iterator):
    '''Abstract base class for database. Requires implementation of make_allele_instance
    Attributes
        db_fasta_path (str)         
        alleles_dict ({})               Dict of {allele_id: AlleleReference instance}. Accessed through direct indexing of allele id (see __getitem__ and __setitem__)
        genes ({})                      Dict of {gene_id: GeneReference instance}
        gap_delimiter (char)            character to denote gaps in aligned sequences                   
        num_landmarks (bool, int)       Flag if allele landmarks have been added, int of number of landmarks if True
    Properties
        alleles                         List of alleles ids
    '''    
    gap_delimiter = '-'
    num_landmarks = None
    allele_distances_path = None

    def __init__(self, db_fasta_path, consensus_path='', ignored_alleles_path='', allele_distances_path=None, gene_clusters_path='', gene_type='IGHV'):
        '''Args
                db_fasta_path (str)             Path to fasta of all alleles. Must be MSA alignment.
                gap_delimiter (str/char)        Char used to indicate gaps in db_fasta_path
                consensus_path (str)            Path to fasta of consensus gainst which to call variants
                gene_clusters_path (str)        Path to TSV of gene's to consider 'same' based on external sequence identity analysis
        '''
        self.gene_type = gene_type.upper()
        # load consensus
        if consensus_path:
            if os.path.exists(consensus_path):
                log.info(f'Using consensus from {consensus_path}')
                consensus = list(SeqIO.parse(consensus_path, 'fasta'))
                if len(consensus) > 1:
                    raise ValueError(f'Consensus fasta {consensus_path} has > 1 sequence')
                self.consensus = str(consensus[0].seq)
            else:
                raise ValueError('\nAllele consensus sequence does not exist at {}'.format(consensus_path))
        else:
            self.consensus = None

        # load gene clusters
        if gene_clusters_path:
            if os.path.exists(gene_clusters_path):
                self.load_gene_clusters(gene_clusters_path)
            else:
                log.warn('\nGene clusters file does not exist at {}'.format(gene_clusters_path))
        else:
            self.gene_clusters = []
            self.gene_in_cluster = {}


        # load ignored alleles
        if ignored_alleles_path:
            if os.path.exists(ignored_alleles_path):
                log.info(f'Using ignored alleles from {ignored_alleles_path}')
                with open(ignored_alleles_path, 'r') as f:
                    self.ignored = set([x.strip() for x in f.readlines()])
            else:
                raise ValueError('\nIgnored alleles file does not exist at {}'.format(ignored_alleles_path))
        else:
            self.ignored = set()

        self.db_fasta_path = db_fasta_path
        alleles_dict = {}
        self.genes = {}
        for record in SeqIO.parse(db_fasta_path, 'fasta'):
            allele = self.make_allele_instance(record.description, record.seq, self.consensus, self.gap_delimiter)
            allele.is_ignored = True if allele.id in self.ignored else False
            alleles_dict[allele.id] = allele
            if allele.gene in self.genes:
                self.genes[allele.gene].add(allele)
            else:
                gene = GeneReference(allele.gene)
                gene.add(allele)
                self.genes[gene.id] = gene
        self.alleles_dict = OrderedDict(sorted([x for x in alleles_dict.items()], key=lambda y:y[0]))

        # set path for allele_distances
        self.allele_distances_path = allele_distances_path if allele_distances_path else os.path.splitext(db_fasta_path)[0].replace('aligned','allele_distance') + '.pickle'

    #
    # Properties
    #
    @property
    def alleles(self):
        return list(self.alleles_dict.values())

    #
    # Magic methods
    #
    def __getitem__(self, key):
        try:
            return self.alleles_dict[key]
        except KeyError:
            raise NotInDatabaseError(f'Allele {key} not in database')
    
    def __setitem__(self, key, value):
        self.alleles_dict[key] = value
    
    def __len__(self):
        return len(self.alleles_dict)

    def __str__(self):
        return f'''Allele database from: {self.db_fasta_path}\n{list(self.alleles_dict.keys())}'''

    def __contains__(self, value):
        return True if value in self.alleles_dict else False

    def __iter__(self):
        yield from self.alleles_dict.keys()
    
    #
    # Abstract methods
    #
    @abstractmethod
    def make_allele_instance(self, *args, **kwargs):
        '''Returns an instance of AlleleReference. 
        Allows for generalization to different fasta formats by implementing subclasses of AlleleReference'''
        raise NotImplementedError

    #
    # Instance methods
    #

    def __next__(self):     # implementing Iterator ABC
        return list(self.alleles_dict.keys()).pop()

    def make_allele_distances(self, allele_distances_path=None):  
        if not allele_distances_path:
            allele_distances_path = self.allele_distances_path
        with open(allele_distances_path, 'rb') as f:
            allele_distances = pickle.load(f)
        if all([x in self.keys() for x in allele_distances]):   # allele_distances dict contains valid allele ids as keys
            self.allele_distances = allele_distances
        
        else: #parsing keys is needed
            new_allele_distances = OrderedDict()

            def parser(key_str):
                result = [x for x in self.keys() if x in key_str]
                if len(result) > 1:
                    raise ValueError(f"Allele distances pickle entry {key_str} matches more than 1 allele id: {str(result)}")
                if len(result) == 0:
                    raise ValueError(f"Allele diostances pickle entry {key_str} does not match any allele id")
                return result[0]
            
            for key, value in allele_distances.items():
                new_allele_distances[parser(key)] = OrderedDict([(parser(k), v) for k, v in value.items()])

            self.allele_distances = new_allele_distances


    def load_gene_clusters(self, path):
        '''Loads TSV of gene clusters (one cluster per line, gene ids tab delimited) and sets attr 
        Sets
            self.gene_clusters                  list of frozenset objects
            self.gene_in_cluster                Dictionary with gene id as key, relevant frozenset cluster as value'''
        clusters = []
        gene_in_cluster = {}
        with open(path, 'r') as f:
            for line in f.readlines():
                cluster = [x.strip() for x in line.split('\t')]
                for gene in cluster:
                    gene_in_cluster[gene] = cluster
                clusters.append(cluster)
        self.gene_clusters = clusters
        self.gene_in_cluster = gene_in_cluster
    
    def build_similar_alleles(self, is_similar=None, variant_filter=None):
        '''Iterates through all pairs of alleles, assigning allele.similar_alleles ([]) using AlleleReference.allele_is_similar'''
        if not variant_filter:
            log.info('Building close alleles by only considering SNPs')
            variant_filter = lambda x: [(pos, var_type) for pos, var_type in x if 'SNP' in var_type]
        if not is_similar:
            threshold = 3
            log.info(f'Allowing at most {threshold} variants for allele pair to be consider as close')
            is_similar = lambda x: True if len(x) <= threshold else False

        for a in self.alleles:
            a.similar_alleles = []
        from math import comb
        num_comb = comb(len(self), 2)
        tenth = round(num_comb/10.0)
        for i, (x, y) in enumerate(combinations(self.alleles, 2)):
            if i % 10000 == 0:
                log.info(f'Proccessed {round(float(i+1)/num_comb*100)}%')
            if x.allele_is_similar(y, is_similar=is_similar, variant_filter=variant_filter):
                x.similar_alleles.append(y)
                y.similar_alleles.append(x)

    def keys(self):
        '''For backwards compatibility with old dictionary allele databases'''
        return self.alleles_dict.keys()
    
    def values(self):
        '''For backwards compatibility with old dictionary allele databases'''
        return self.alleles_dict.values()
    
    def iteritems(self):
        return list(self.alleles_dict.items())
    def items(self):
        return list(self.alleles_dict.items())
    
    def make_landmarks(self, num_landmarks, read_length, sampled_depth, sampled_variance, sampled_edge_variance, minimum_coding_bases, num_landmark_groups=1, max_copies=4):
        '''
        Args
            sampled_variance                Single copy sampled variance
            sampled_depth                  Single copy sampled sequencing depth
            sampled_edge_variance          Sample from TTN using bam_filter_classes    [mean pos 1 sampled variance,..., mean pos minimum_coding_bases sampled variance]'''
        log.info(f"Making allele landmark positions (n={num_landmarks}, groups={num_landmark_groups})")
        
        self.num_landmarks = num_landmarks
        self.num_landmark_groups = num_landmark_groups
        self.read_length = read_length
        self.sampled_depth = sampled_depth
        self.sampled_variance = sampled_variance
        self.minimum_coding_bases = minimum_coding_bases
        self.depth_coefficient = 1/(self.read_length/float(self.sampled_depth))
        self.max_copies = max_copies
        self.sampled_edge_variance = sampled_edge_variance

        for allele in self.alleles_dict.values():
            allele.landmark_groups, allele.landmarks = self.get_allele_landmarks(allele)        


    def get_allele_landmarks(self, allele):
        '''Chooses landmark positions, gets expected depth and stdev from calculate_landmark_expection, assigns landmarks to groups, returns groups (list of lists) and landmarks (list)'''
        num_landmarks = self.num_landmarks+1
        a_len = len(allele)
        
        # if a_len < self.minimum_coding_bases:   # allele is too short to have valid mappings
        #     return None, None
        
        step = round(a_len/float(num_landmarks))
        landmarks = [x for x in range(step, a_len, step)][:self.num_landmarks]
                
        landmark_groups = [[] for i in range(self.num_landmark_groups)]
        for i, pos in enumerate(landmarks):          
            l = self.calculate_landmark_expection(allele, pos)
            try:
                landmark_groups[i%self.num_landmark_groups].append(l)
            except IndexError:
                log.error((i, len(landmarks), self.num_landmark_groups, i%self.num_landmark_groups, len(landmark_groups)))
                raise
            landmarks[i] = l

        return landmark_groups, landmarks
    
    def calculate_landmark_expection(self, allele, pos):
        '''Returns a tuple representing an allele landmark: (position, expected_depth, [0, 1 copy expected stddev, 2 copy stdev,..., self.max_copies expected stdev])'''
        allele.last_start = len(allele)-self.minimum_coding_bases
        flanking_length = self.read_length - self.minimum_coding_bases

        relative_pos = min(pos, len(allele)-pos)

        num_start_positions = min(self.read_length, flanking_length+pos+1) - max(0, pos-allele.last_start)
        single_copy_trial_probability = (float(self.sampled_depth))/self.read_length
        
        expected_depth = num_start_positions*single_copy_trial_probability

        stdevs = [0]
        for copy_num in range(1, self.max_copies+1):
            if relative_pos < self.minimum_coding_bases:
                stdevs.append(sqrt(copy_num*self.sampled_edge_variance[relative_pos]))
            else:
                stdevs.append(sqrt(self.sampled_variance*copy_num))
        
        return (pos, int(round(expected_depth)), stdevs)
    
    def get_expected_depth(self, allele, pos):
        '''Some duplicate code from calculate_landmark_expection but used in analysis functions'''
        allele.last_start = len(allele)-self.minimum_coding_bases
        flanking_length = self.read_length - self.minimum_coding_bases

        num_start_positions = min(self.read_length, flanking_length+pos+1) - max(0, pos-allele.last_start)
        single_copy_trial_probability = (float(self.sampled_depth))/self.read_length
        
        expected_depth = num_start_positions*single_copy_trial_probability
        return expected_depth

    @property
    def novel_allele_class(self, *args, **kwargs):
        return NovelAlleleReference(*args, **kwargs)
    
    def add_novel_allele(self, allele_id, *args, **kwargs):
        '''Creates a new allele object instance as defined by novel_allele_class property, adds to appropriate attributes'''
        allele = self.novel_allele_class(allele_id, *args, **kwargs)
        
        self.alleles_dict[allele.id] = allele
        if allele.gene in self.genes:
            self.genes[allele.gene].add(allele)
        else:
            gene = GeneReference(allele.gene)
            gene.add(allele)
            self.genes[gene.id] = gene


    def add_allele_clusters(self, cluster_path):
        '''Takes a csv of allele clusters (once cluster per line) and adds the list of alleles in the cluster as a attribute .cluster to
        every allele object in the database'''
        from .allele_clusters import Cluster
        self.clusters = []
        self.allele_to_cluster = defaultdict(lambda: None)
        for a in self.alleles: a.cluster = None
        with open(cluster_path, 'r') as f:
            for i, line in enumerate(f.readlines()):
                alleles = [self[a.strip()] for a in line.split('\t')]
                cluster = Cluster(i, alleles)
                for allele in alleles:
                    allele.cluster = cluster
                    self.allele_to_cluster[allele.id] = cluster
                self.clusters.append(cluster)

    def get_wildtype_target(self, gene: str):
        ''' In many cases this returns the *01 allele, except with genes that are present in self.gene_in_cluster
        in which case it returns the first allele of the first gene of the alphabetical sort of the cluster '''
        if gene not in self.gene_in_cluster:
            if gene+'*01' in self:
                wildtype_target=gene+'*01'
            else:
                wildtype_target = sorted(self.genes[gene].alleles, key=lambda x: x.id)[0].id
        else:
            cluster = self.gene_in_cluster[gene]
            wildtype_target = sorted(self.genes[cluster[0]].alleles, key=lambda x: x.id)[0].id # should give *01 in most cases
        return wildtype_target


    def print_alignment(self, a, b, bp_per_line=130):
        output = []
        a_seq = self[a].aligned_seq
        b_seq = self[b].aligned_seq
        alignment = self[a].alignment(self[b])
        for i in range(0, len(a_seq), bp_per_line):
            output.append(a_seq[i:i+bp_per_line]+'\t'+a)
            output.append(alignment[i:i+bp_per_line])
            output.append(b_seq[i:i+bp_per_line]+'\t'+b)
            output.append(f"\n{'- '*int(bp_per_line/2)}\n")
        print('\n'.join(output))
    
    def add_novel_allele(self, allele_id, aligned_seq, functional):
        instance = NovelAlleleReference(allele_id, aligned_seq, functional, consensus_seq=self.consensus, gap_delimiter=self.gap_delimiter)
        self.alleles.append(instance)
        self.alleles_dict[instance.id] = instance

class ImgtAlleleDatabase(AlleleDatabase):
    '''Child class of AlleleDatabase that uses ImgtAlleleReference class for alleles
    Attributes
        gap_delimiter        Set as ImgtAlleleReference in __init__
    '''
    gap_delimiter = '-'
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def make_allele_instance(self, *args, **kwargs):
        return ImgtAlleleReference(*args, **kwargs)


class ImgtNovelAlleleDatabase(AlleleDatabase):
    '''Child class of AlleleDatabase that uses ImgtAlleleReference class and NovelAlleleReference for alleles
    Attributes
        gap_delimiter        Set as ImgtAlleleReference in __init__
    '''
    gap_delimiter = '-'
    def __init__(self, *args, flanking_json=None, minimum_mapped_prefix_bases=100, **kwargs):
        if flanking_json:
            with open(flanking_json, 'r') as f:
                self.flanking_data = json.load(f)
        else:
            self.flanking_data = None

        super().__init__(*args, **kwargs)
    
    def make_allele_instance(self, description, seq, consensus, gap_delimiter):
        if 'Homo_sapiens' in description: # allele is IMGT
            return ImgtAlleleReference(description, seq, consensus, gap_delimiter, self.flanking_data)
        elif 'Novel' in description:   # is oscar novel allele
            return OscarNovelAlleleReference(description, seq, consensus, gap_delimiter)
        else:
            raise NotImplementedError(f"Allele type for {description} not implemented")

    def calculate_landmark_expection(self, allele, pos):
        '''Returns a tuple representing an allele landmark: (position, expected_depth, [0, 1 copy expected stddev, 2 copy stdev,..., self.max_copies expected stdev])'''

        if not allele.has_flanking:
            return super().calculate_landmark_expection(allele, pos)
        
        # check if relevant flanking is long enough to not scale depth
        start_closest = True if pos < (len(allele)/2) else False
        if (start_closest and  len(allele.upstream_flanking)+pos >= self.read_length) or (not start_closest and len(allele.downstream_flanking)+(len(allele)-pos) >= self.read_length):
            return (pos, round(self.sampled_depth), [0]+[sqrt(self.sampled_variance * copy_num) for copy_num in range(1, self.max_copies+1)])
        

        # landmark pos close enough to beginning/end that it needs to be scaled
        relative_pos = min(pos, len(allele)-pos)
        flanking_len = len(allele.upstream_flanking) if start_closest else len(allele.downstream_flanking)

        expected_depth_scaling_factor = min((self.read_length-self.minimum_coding_bases+flanking_len+relative_pos)/self.read_length, 1)
        
        expected_depth = expected_depth_scaling_factor * self.sampled_depth
        single_copy_variance = expected_depth_scaling_factor * self.sampled_variance     # scale single copy variance to get landmark single copy variance

        stdevs = [0] + [sqrt(copy_num * single_copy_variance) for copy_num in range(1, self.max_copies+1)]       # scale landmark single copy variance by copy number
        
                
        return (pos, round(expected_depth), stdevs)
    

    @staticmethod
    def get_edit_distance(a1, a2, ignore_start_gaps: bool=False, ignore_end_gaps: bool=False) -> int:
        """Calculates edit distance using AlleleRferences.aligned_seq
        Args:
            ignore_start_gaps:          does not include prefix gaps in edit distance output
            ignore_end_gaps:            does not include suffix gaps in edit distance output         
        """

        first_m_flag = False
        leading_gaps = 0
        gaps_since_last_m = 0
        edit_distance = 0
        for i, j in zip(a1, a2):
            if i != j:
                edit_distance += 1
                if (i == '-' or j == '-'):
                    if not first_m_flag:
                        leading_gaps += 1
                    gaps_since_last_m += 1
                else:
                    gaps_since_last_m = 0
                    first_m_flag = True
                
            else:
                gaps_since_last_m = 0
                first_m_flag = True

        return edit_distance - (ignore_start_gaps * leading_gaps) - (ignore_end_gaps * gaps_since_last_m)

    


#
# Allele and Gene classes
#

class AlleleReference(ABC):
    ''' Abstract base class for allele references. Requires implementation of set_attr
    Attributes
        id (str)                Allele name
        aligned_seq (str)       Aligned allele sequence with gap_delimiter as gaps
        functional (str)        'F'=functional, 'P'=pseudogene, 'ORF'=open reading frame
        is_ignored (Bool)       Boolean to indicate if the allele is an 'ignored' allele. None if not specified.
        gap_delimiter (char)    Char used to indicate gaps in db_fasta_path
        variants (set)          Set of variants that define the allele. Called against the provided consensus sequence. None of consensus not provided.
                                !!! Variants are relative to the *unaligned positions* on the *consensus* sequence !!!
    Properties
        is_functional           Returns true if self.functional == 'F' else False
    '''

    def __init__(self, allele_data, aligned_seq, consensus_seq=None, gap_delimiter='-', *args, **kwargs):
        self.aligned_seq = str(aligned_seq)
        self.seq = self.aligned_seq.replace(gap_delimiter,'')
        self.set_attr(allele_data, *args, **kwargs)
        self.gap_delimiter = gap_delimiter
        if consensus_seq:
            self.set_variants(consensus_seq)
        else:
            self.variants = None

    #
    # Properties
    #
    @property
    def is_functional(self):
        return True if self.functional == 'F' else False

    @property
    def reverse_complement(self):
        return str(Seq(self.seq).reverse_complement())
    
    @property
    def has_flanking(self):
        try:
            return self._has_flanking
        except AttributeError:
            return False
    
    #
    # Magic methods
    #
    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return f'''{self.id}
Length:     {len(self)}
Functional: {self.functional}'''

    def __getitem__(self, key):
        '''Essentially a getter to ensure backwards compatibility when the allele db was a dict of dicts'''
        if key == 'variants':
            return self.variants
        if key == 'length':
            return len(self)
        if key == 'functional':
            return self.functional
        if key == 'seq':                # this is because 'seq' in the old db was the aligned seq
            return self.aligned_seq
        else:
            raise KeyError(f'Provided key "{key}" not a attribute or not implemented')

    #
    # Abstract methods
    #
    @abstractmethod
    def set_attr(self, data, *args, **kwargs):
        raise NotImplementedError

    #
    # Instance methods
    #
    def allele_is_similar(self, other_allele, 
                                is_similar=lambda x: True if len(x) < 5 else False,
                                variant_filter=None):
        '''Builds variants with other_allele (AlleleReference instance), returns True if the filtered variants (using variant_filter) matches is_similar'''
        if not variant_filter:
            variant_filter = lambda x: len(self.get_snps(other_allele))
        variants = variant_filter(self.get_variants(self.aligned_seq.replace(self.gap_delimiter, '.'), other_allele.aligned_seq.replace(other_allele.gap_delimiter, '.'))) 
        return is_similar(variants)
    
    def get_snps(self, other_allele):
        return [(pos, var_type) for pos, var_type in 
                self.get_variants(other_allele.aligned_seq.replace(other_allele.gap_delimiter, '.'), self.aligned_seq.replace(self.gap_delimiter, '.'))
                if 'SNP' in var_type]

    def get_all_variants(self, other_allele):
        return self.get_variants(other_allele.aligned_seq.replace(other_allele.gap_delimiter, '.'), self.aligned_seq.replace(self.gap_delimiter, '.'))

    def num_snps(self, other_allele):
        return len(self.get_snps(other_allele))

    def set_variants(self, consensus_seq):
        '''Driver to use get_variants, remove_periods'''
        try:
            self.variants = self.get_variants(self.seq.replace(self.gap_delimiter, '.'), consensus_seq.replace(self.gap_delimiter, '.'))
        except ValueError:
            raise ValueError(f'Errors while building variants for {self.id} (see debug log)')

    def mismatches(self, other_allele):
        '''Provided an instance of AlelleReference returns (number of mismatches, string where ' ' is match, '*' is mismatch)'''
        distance = 0
        alignment = ''
        for i, _ in enumerate(self.aligned_seq):
            if self.aligned_seq[i] != other_allele.aligned_seq[i]:
                distance += 1
                alignment = alignment+'*'
            else:
                alignment = alignment+' '
        return distance, alignment
    
    def alignment(self, other_allele):
        return self.mismatches(other_allele)[1]
    

    #
    # Static methods - these are old implentations, sorry for the mess
    #

    @staticmethod
    def get_variants(allele_seq, consensus_seq):
        ## Calls variants for allele_seq relative to consensus_seq
        ## Input:	Input sequences must be aligned with '.' as the inbdel character
            
        variants = []
        msg = []

        for index, pos in enumerate(consensus_seq):         
            try:
                same_value = (allele_seq[index] == pos)
                # log.debug('Index: {} seq: {} len: {} consensus: {} len: {}'.format(index, pos, len(allele_seq), consensus_seq[index], len(consensus_seq)))

            except IndexError as e:
                log.debug("Allele seq is shorter than consensus, not adding 3' deletion")
                break

            if same_value:
                continue
            elif allele_seq[index] == '.':  # allele has deletion relative to consensus
                if (variants and 'DEL' in variants[-1]['op'] and                                    # check if previous variant was deletion and check that all nucl since start of deletion were .
                            all(x == '.' for x in allele_seq[variants[-1]['pos']:index])):  
                    variants[-1]['op'] = variants[-1]['op'] + consensus_seq[index]          
                else:
                    variants.append({'pos': index, 'op': 'DEL.{}'.format(pos)})
            elif consensus_seq[index] == '.':                                                       #allele has an insertion relativec to consensus
                if (variants and 'INS' in variants[-1]['op'] and                                    # check if previous variant was insertions and check that all nucl since start of deletion were .
                            all(x == '.' for x in consensus_seq[variants[-1]['pos']:index])):
                    variants[-1]['op'] = variants[-1]['op'] + allele_seq[index]                           
                else:
                    variants.append({'pos': index, 'op': 'INS.{}'.format(allele_seq[index])})
            else:
                variants.append({'pos': index, 'op': 'SNP.{}{}'.format(pos, allele_seq[index])})
        
        # ## Remove varaints that are due to truncated 3' or 5' in allele_seq or consensus_seq
        # variants = [x for x in variants if not ('DEL' in x['op'] and (x['pos'] == 0 or x['pos']+len(x['op'].split('.')[1]) >= len(allele_seq)-1))]
        # variants = [x for x in variants if not ('INS' in x['op'] and ((x['pos'] == 0 and consensus_seq[0] == '.') or (x['pos']+len(x['op'].split('.')[1]) >= len(allele_seq)-1 and consensus_seq[-1] == '.')))]
        # log.debug(variants)
        consensus_seq_np, forward, backward = AlleleReference.remove_periods(consensus_seq)
        converted_variants = AlleleReference.adjust_variant_positions(variants, forward, consensus_seq_np, allele_seq.replace('.',''))
        
        return 	set([(x['pos'], x['op']) for x in converted_variants])

    @staticmethod
    def remove_periods(seq):
        # returns (seq with '.' removed, array of len(result) where each pos maps to corresponding pos in original seq)
        forward_map = []
        backward_map = []
        result = []
        seq = seq.replace('-', '.')
        for index, char in enumerate(seq):
            if char != '.':
                result.append(char)
                backward_map.append(index)
                forward_map.append(len(result)-1)
            else:
                forward_map.append(0 if not forward_map else (forward_map[-1]+1 if seq[index-1] != '.' else forward_map[-1]))
        return (''.join(result), forward_map, backward_map)
    
    @staticmethod
    def adjust_variant_positions(variants, forward_mapping, consensus_seq_np, allele_seq_np):
        def test_variants(variants, consensus_seq, allele_seq, not_converted_variants):
            ## Tests variant accurancy by applying variants to consensus sequence and comparing to allele sequence
            ## Input: 	variants = [{'pos': int, 'op': VTYPE.NN}, ....]
            ##			consensus_seq, allele_seq must have periods removed
            ## Output: [] if no errors, o/w [(log.messagetype, 'Message text {}')]

            result = []
            for var in variants:
                if 'SNP' in var['op'] and consensus_seq[var['pos']] != var['op'].split('.')[1][0]:
                    result.append((log.error, '{} at position {} does not correspond with consensus seq'.format(var['op'], var['pos'])))
            translate = ''
            var_index = 0
            index = 0
            while index < len(consensus_seq):
                if var_index < len(variants) and index == variants[var_index]['pos']:
                    var = variants[var_index]
                    if 'SNP' in var['op']:
                        translate = translate + var['op'].split('.')[1][1]
                    elif 'DEL' in var['op']:
                        index += len(var['op'].split('.')[1]) - 1
                    elif 'INS' in var['op']:
                        translate = translate + var['op'].split('.')[1] 
                        index -= 1
                    if translate != allele_seq[:len(translate)]:
                        result.append((log.error, 'Variant {} is incorrect\nBefore conversion: {}\nAllele seq\n{}\n{}\nConverted seq'.format(var, not_converted_variants[var_index], allele_seq[:len(translate)], translate)))
                        break
                    var_index += 1
                else:
                    translate += consensus_seq[index]
                index += 1
            return result

        ## For converting variant position from aligned to un-aligned (ie w/o periods) positions
        ## Input:	forward_mapping is forward_map output of remove_periods()
        converted_variants = []
        for v in variants:
            # new_pos = (forward_mapping[v['pos']] if 'INS' not in v['op'] else forward_mapping[v['pos']]+1)
            new_pos = forward_mapping[v['pos']]
            converted_variants.append({'pos': new_pos, 'op': v['op']})

        msg = test_variants(converted_variants, consensus_seq_np, allele_seq_np, variants)
        if msg:
            for f, message in msg:
                print(message)
                f(message)
            raise ValueError('Error in variant calling - see log messages')

        return converted_variants


class ImgtAlleleReference(AlleleReference):
    ''' Child class implementation of AlleleReference for IMGT-style fasta allele format
    Attributes
        accession (str)                 IMGT/LIGM-DB accession number(s)
        accession_start (int)           start position in the IMGT/LIGM-DB accession number(s)
        accession_end = (int)           end position in the IMGT/LIGM-DB accession number(s)
        accession_len (int)             number of nucleotides in the IMGT/LIGM-DB accession number(s)
        coding_start (int)              codon start, or None for non coding labels
        upstream_bases (int)            +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
        downstream_bases (int)          +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
        is_upstream_partial = False     Sequence is partial in 5' 
        is_downstream_partial = False   Sequence is partial in 3'
        original_imgt_description       Source IMGT-formatted allele description
    Properties
        imgt_description                Generates the IMGT-formatted allele description using current value of attributes
    '''
    @property
    def imgt_description(self):
        partial = '_'
        if self.is_upstream_partial or self.is_downstream_partial:
            partial = 'partial_in_'
            if all([self.is_upstream_partial, self.is_downstream_partial]):
                partial = partial + "3'_and_in_5'"
            else:
                partial = partial + ("3'" if self.is_downstream_partial else "5'")

        return '|'.join([str(x) for x in [self.accession, 
                        self.id,
                         'Homo_sapiens', 
                         self.functional, 
                         'V-REGION',
                         f'{self.accession_start}..{self.accession_end}',
                         self.accession_len,
                         self.coding_start if self.coding_start else 'NR',
                         self.upstream_bases if self.upstream_bases else '_',
                         self.downstream_bases if self.downstream_bases else '_',
                         '_',
                         '_',
                         '_',
                         partial,
                         'rev-compl' if self.is_reversed else '_']])+'|'


    def set_attr(self, data, flanking_data=None):
        self.original_imgt_description = data
        data = data.split('|')
        try:
            self.id = data[1]
            self.gene = self.id.split('*')[0]
            
            self.accession = data[0]
            self.id = data[1]
            self.functional_raw = data[3]
            self.functional = self.functional_raw[1] if self.functional_raw[0] in {'[', '('} else self.functional_raw
            accession_interval = data[5]
            self.accession_start, self.accession_end = [int(x) for x in data[5].split('..')] if '..' in accession_interval else (None, None)
            self.accession_len = int(data[6].split(' ')[0].split('_')[0])
            try:
                self.coding_start = int(data[7])
            except ValueError:
                if data[7] == 'NR': 
                    self.coding_start = None
                else:
                    raise ValueError(f'Value {data[7]} for coding start is invalid\nAllele entry:\n{"|".join(data)}')
            self.upstream_bases = None if not data[8] or data[8] in [' ', '_'] else int(data[8])
            self.downstream_bases = None if not data[9] or data[9] in [' ', '_'] else int(data[9])
            
            partial = data[13]
            self.is_upstream_partial = True if "5'" in partial else False
            self.is_downstream_partial = True if "3'" in partial else False

            self.is_reversed = True if 'rev-compl' in data[14] else False

            # Add flanking sequence data
            if flanking_data and (self.id in flanking_data) and (flanking_data[self.id]['start'] > 0 or flanking_data[self.id]['end'] > len(flanking_data[self.id]['sequence'])):
                self._has_flanking = True
                self.flanking_data = flanking_data[self.id]
                if self.flanking_data['sequence'][self.flanking_data['start']:self.flanking_data['end']] != self.seq:
                    raise ValueError(f"Sequence in flanking sequence json for allele {self.id} does not match allele sequence")
                self.upstream_flanking = self.flanking_data['sequence'][:self.flanking_data['start']]
                self.downstream_flanking = self.flanking_data['sequence'][self.flanking_data['end']:]
                self.seq_with_flanking = self.flanking_data['sequence']
                self.coding_start = self.flanking_data['start']
                self.coding_end = self.flanking_data['end']
            else: # flanked by Ns
                self._has_flanking = True
                self.upstream_flanking = 'n'*100
                self.downstream_flanking = 'n'*100
                self.seq_with_flanking = ('n'*100) + self.seq + ('n'*100)
                self.coding_start = 100
                self.coding_end = len(self.seq)-100



        except (ValueError, IndexError) as e:
            raise ValueError(f'Allele data not valid {data}').with_traceback(e.__traceback__)
    
class NovelAlleleReference(AlleleReference):
    ''' Child class implementation of AlleleReference for novel alleles
    Attributes
    Properties
    '''

    def __init__(self, allele_id, aligned_seq, functional, consensus_seq, gap_delimiter='-'):
        '''
        Args
            data (dict):        key-value pairs of attributes to set'''
        self.id = allele_id
        self.aligned_seq = str(aligned_seq)
        self.seq = self.aligned_seq.replace(gap_delimiter,'')
        self.gap_delimiter = gap_delimiter
        if consensus_seq:
            self.set_variants(consensus_seq)
        else:
            self.variants = None
        self.functional = functional

        self._has_flanking = True
        self.upstream_flanking = 'n'*100
        self.downstream_flanking = 'n'*100
        self.seq_with_flanking = ('n'*100) + self.seq + ('n'*100)
        self.coding_start = 100
        self.coding_end = len(self.seq)-100



    def set_attr(self, seq, functional_type):
        pass

class OscarNovelAlleleReference(AlleleReference):
    ''' Child class implementation of AlleleReference for novel alleles
    Attributes
    Properties
    '''

    def set_attr(self, data):
        self.id = data
        self.gene = self.id.split('*')[0]
        self.functional = 'F' if "or" not in self.id else "P"

        self._has_flanking = True
        self.upstream_flanking = 'n'*100
        self.downstream_flanking = 'n'*100
        self.seq_with_flanking = ('n'*100) + self.seq + ('n'*100)
        self.coding_start = 100
        self.coding_end = len(self.seq)-100



class GeneReference(object):
    ''' Class to group / keep track of genes and their associated alleles
    Attributes
        id (str)            Gene id
        alleles ([])        list of AlleleReference objects
    '''

    def __init__(self, id, alleles=[]):
        self.id = id
        self.alleles = []
        if isinstance(alleles, AlleleReference):
            self.alleles.add(alleles)
        elif isinstance(alleles, list):
            for a in alleles:
                self.alleles.add(a)
        else:
            raise ValueError('alleles argument requires instance or list of AlleleReference objects')

    #
    # Magic methods
    #
    def __eq__(self, other):
        return True if self.id == str(other) else False
    
    def __str__(self):
        return self.id
    
    def __len__(self):
        return len(self.alleles)

    #
    # Attributes
    # 

    @property
    def is_functional(self):
        return any([a.is_functional for a in self.alleles])

    #
    # Instance methods
    #
    def add(self, allele):
        '''Adds allele to self.alleles
        Args
            allele              Instance of AlleleReference object
        '''
        if not isinstance(allele, AlleleReference):
            raise ValueError('Allele must be instance of AlleleReference')    
        self.alleles.append(allele)
    
class NotInDatabaseError(Exception):
    """Raise instead of key error when provided allele is not in database"""
    pass
