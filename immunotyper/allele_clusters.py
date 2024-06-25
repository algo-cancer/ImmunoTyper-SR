from typing import List
from statistics import mean
from collections import Counter, defaultdict

def make_cluster_consensus(cluster: List):
    cluster_seqs = [allele_db[a].aligned_seq for a in cluster]
    consensus = [Counter(pos).most_common()[0][0] for pos in zip(*cluster_seqs)]
    return ''.join(consensus)


class Cluster(object):
    
    def __init__(self, name, alleles, non_mapping_pos):
        self.id = name
        self.alleles = alleles
        self.non_mapping_pos = non_mapping_pos
        self.alleles_dict = dict([(a.id, a) for a in alleles])
        self.consensus_aligned = self.make_cluster_consensus()
        self.consensus, self.align_to_consensus, _ = self.alleles[0].remove_periods(self.consensus_aligned)
        self.consensus_depth = [0 for i in range(len(self.consensus))]
        self.reads = []
    
    def convert_allele_coordinates_to_consensus(self, allele_positions, allele):
        _, _, allele_to_align = self.alleles[0].remove_periods(self.alleles_dict[allele].aligned_seq)
        aligned_positions = [allele_to_align[x] for x in allele_positions]
        consensus_positions = [self.align_to_consensus[x] for x in aligned_positions]
        return consensus_positions
    
    def add_read_depth(self, start, end):
        if start < 0 or end >= len(self.consensus):
            raise ValueError(f"Invalid coordinates for {(start, end)} (len consensus = {len(self.consensus)})")
        
        for i in range(start, end):
            self.consensus_depth[i] += 1

    def get_expected_depth(self, pos, sampled_depth, minimum_coding_bases=50, read_length=150):
        '''Some duplicate code from calculate_landmark_expection but used in analysis functions'''
        last_start = len(self.consensus)-minimum_coding_bases
        flanking_length = read_length - minimum_coding_bases

        num_start_positions = min(read_length, flanking_length+pos+1) - max(0, pos-last_start)
        single_copy_trial_probability = (float(sampled_depth))/read_length

        expected_depth = num_start_positions*single_copy_trial_probability
        return expected_depth

    def set_cnv_call(self, single_copy_depth):
        self.cnv_call = round(mean(self.consensus_depth[pos]/self.get_expected_depth(pos, single_copy_depth) for pos in self.non_mapping_pos))

    def make_cluster_consensus(self):
        cluster_seqs = [a.aligned_seq for a in self.alleles]
        consensus = [Counter(pos).most_common()[0][0] for pos in zip(*cluster_seqs)]
        return ''.join(consensus)


def calculate_read_depths(allele_to_cluster, reads):
    for read in reads:
        already_added = set()
        for allele, mapping in read.mappings_dict.items():
            c = allele_to_cluster[allele]
            if c and c.id not in already_added:
                already_added.add(c.id)
                s, e = c.convert_allele_coordinates_to_consensus((read.get_start(allele), read.get_end(allele)-1), allele)
                c.add_read_depth(s, e-1)
                c.reads.append(read)
