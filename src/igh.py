import sys
import create_allele_db, clustering, msa_generators, cluster_merging, candidate_calling, cluster_breaker
from common import log, SeqRecord, initialize_logger
from Bio import SeqIO
from minimap_wrapper import MinimapWrapper
from cluster_class import SubCluster



def run(reads_path, expected_coverage, blasr_src='blasr', minimap='minimap', spoa='../tools/spoa/build/bin/spoa', poa='../tools/poaV2/poa', dsf_src=None, max_copy=4, cache_prefix='.', no_cov_estimation=False, skip_extraction=False, debug_log_path=None, output_solver_log=False):
	initialize_logger(debug_log_path)
	allele_db = create_allele_db.iterate_fasta()

	# load reads
	reads = fix_rev_comp(reads_path, skip_extraction=skip_extraction, mapper=MinimapWrapper(minimap_src=minimap))
	for index, read in enumerate(reads):
		read.name = read.id
		read.id = str(index)
	log.info('Using {} reads'.format(len(reads)))

	# generate mapping-based clusters
	log.info('-----\nGenerating mapping-based assignments...')
	mapping_cluster_tool = clustering.MappingClustering(expected_coverage)
	mapping_clusters, ambiguous_reads, coverage = mapping_cluster_tool.cluster(reads, no_cov_estimation=no_cov_estimation)
	log.info('Generated {} clusters using mapping\n{} reads identified as ambiguous'.format(len(mapping_clusters), len(ambiguous_reads)))

	# generate clusters for ambiguous reads
	log.info('-----\nClustering ambiguous reads...')
	min_cluster_size = int(0.75*coverage)	# for dsf clustering parameter
	clustering_tool =  clustering.DSF(dsf_src=dsf_src,
										distance_calculator=clustering.MinimapDistanceCalculator(minimap=MinimapWrapper(minimap_src=minimap, params='-cx ava-pb -k14 -w3'), 
																								filter_function=lambda x: True if x <0.3 else False),
										params='--min-size {} --min-fillin 0.95 -t 1'.format(min_cluster_size))
	clusters = clustering_tool.cluster(ambiguous_reads)
	cluster_is_invalid = lambda x: True if len(x) < min_cluster_size else False
	valid_clusters = [c for c in clusters if not cluster_is_invalid(c)]
	log.info('Generated {} valid clusters and {} small cluster using sequence similarity clustering from ambiguous reads'.format(len(valid_clusters), len(clusters)-len(valid_clusters)))

	# merge resultant invalid / single-read clusters with valid clusters
	log.info('-----\nMerging invalid clusters into valid clusters...')
	consensus_builder = msa_generators.Spoa(src=spoa)
	for c in valid_clusters:
		c.consensus = consensus_builder.generate_consensus(c)
	remapper = cluster_merging.ReMapper(to_remap=cluster_is_invalid, mapper=MinimapWrapper(minimap_src=minimap, params="-cx map-pb -k10 -w3 -N5"))
	clusters = remapper.merge_clusters(clusters)

	for c in clusters:
		c.consensus = consensus_builder.generate_consensus(c)

	# break clusters into subclusters using LP model
	log.info('-----\nBreaking clusters...')
	breaker = cluster_breaker.IlpBreaker(allele_db=allele_db, 
										read_depth=coverage,
										max_copy=max_copy,
										poa=poa,
										poa_cache=lambda x: '../cache/{}/{}.poa'.format(cache_prefix, x.id), 
										gurobi_log_path=(None if not (debug_log_path and output_solver_log) else debug_log_path+'.solver.log'))
	clusters = breaker.break_clusters(clusters)
	log.info('{} clusters after breaking'.format(len(clusters)))

	# merge unbroken (ie infeasible or no candidate clusters) into successful subclusters
	number_tries = 0 # loop counter to break after threshold
	cluster_is_invalid = lambda x: True if not isinstance(x, SubCluster) else False
	allele_db = create_allele_db.iterate_fasta()
	remapper = cluster_merging.RemapToParent(to_remap=cluster_is_invalid, mapper=MinimapWrapper(minimap_src=minimap, params="-cx map-pb -k10 -w3 -N5"))
	breaker.to_break = cluster_is_invalid

	while any([cluster_is_invalid(c) for c in clusters]) or number_tries > 2:
		if number_tries: log.info('-----\nAttempt # {} at cluster merging and rebreaking\n------'.format(number_tries))
		log.info('Merging {} clusters that could not be broken'.format(len([c for c in clusters if cluster_is_invalid(c)])))
		
		for c in clusters:
			if isinstance(c, SubCluster):		# therefor successful subcluster
				c.consensus = consensus_builder.generate_consensus(c)
		clusters = remapper.merge_clusters(clusters)

		# re-break newly merged clusters
		clusters_to_break = [x for x in clusters if cluster_is_invalid(x)]
		log.info('-----\nRebreaking {} newly merged clusters'.format(len(clusters_to_break)))
		for c in clusters_to_break:
			c.consensus = consensus_builder.generate_consensus(c)
		breaker.poa_cache = lambda x: '../cache/{}/merged/{}.poa'.format(cache_prefix, x.id)
		clusters = breaker.break_clusters(clusters, maxcopy=max_copy)
		log.info('{} clusters after breaking'.format(len(clusters)))

		number_tries += 1

	# finally call alleles
	log.info('Calling alleles for clusters')

	for c in mapping_clusters:
		c.consensus = consensus_builder.generate_consensus(c)

	allele_caller = candidate_calling.BlasrCaller(blasr_src=blasr_src, consensus_builder=consensus_builder)
	map(allele_caller.call_cluster, clusters)

	return mapping_clusters+clusters


def main():
	reads_path = sys.argv[1].strip()
	expected_coverage = int(sys.argv[2].strip())
	dsf_src = sys.argv[3].strip() if len(sys.argv) > 3 else None

	run(reads_path, expected_coverage, dsf_src)


def fix_rev_comp(reads_path, 
				 mapper=MinimapWrapper(),
				 allele_db='../database/V-QUEST-reference-allele-db+no-period-references.clustalw.no-gaps.fasta',
				 num_mappings_to_save=5, skip_extraction=False):
	
	mapper.params = '-cx map-pb -k10 -w3 -N{}'.format(num_mappings_to_save-1)
	try:
		mappings = mapper.run(reads_path, allele_db)
		reads = dict([(x.id, x) for x in SeqIO.parse(reads_path, 'fasta')])
	except IOError, ValueError:
		log.error('Reads file does not exist or is invalid')
		raise ValueError

	unmapped = set(reads.keys()).difference([x.qName for x in mappings])
	if len(unmapped) > 0:
		log.warn('{} reads had no allele mapping, they will be removed\n'.format(len(unmapped), '\n'.join(list(unmapped))))
		log.debug('Read ids to be removed:\n'+'\n'.join(unmapped))
		reads = dict([(k, v) for k, v in reads.iteritems() if k not in unmapped])
	log.info('Loaded {} reads'.format(len(reads)))
	
	if not skip_extraction: reads, mappings = get_read_segments(reads, mappings)

	modified_seqs = []
	for m in mappings:
		try:
			reads[m.qName].mapping.append(m)
		except AttributeError:
			reads[m.qName].mapping = [m]
	for r in reads.values():
		r.mapping = sorted(r.mapping, key=lambda x: x.tErrors())
		if r.mapping[0].strand == '-':
			r.seq = r.seq.reverse_complement()
			modified_seqs.append(r)
	log.info('{} reads were reverse complemented'.format(len(modified_seqs)))
	
	
	return reads.values()

def get_read_segments(reads, mappings, flanking=1000):
	"""Extract segments of reads containing IGHV genes plus flanking sequence. Multiple segments from the same read have 

	Args:
		reads: 					Dictionary of original reads
		mappings:				List of minimap_wrapper.PAFMapping objects encoding mapping of reads to V gene database
	Returns:
		x: 						Dict of (id, SeqIO.SeqRecord-like read object) containing read segments
		mappings: 				Input mappings but with the qNames updated to the read segment ids
	"""
	log.info('Extracting read segments with {} flanking bp'.format(flanking))
	extracted_reads = {}
	read_mappings = []
	curr_read = mappings[0].qName
	clipped_segments = []
	read_segments = []
	no_valid_segments = []
	new_mappings = []
				  
	is_valid = lambda x: True if x.tEnd-x.tStart > 0.9*x.tLength else False		# mapping is valid if covers >90% of IGHV gene sequence
	for m in mappings:
		if m.qName != curr_read: # finished with mappings from current read
			if read_mappings: # save segments
				curr_read_segments = []
				read = reads[curr_read]
				extracted_reads[curr_read] = []
				for i, (start, end, mapping) in enumerate(read_mappings):
					if (start < flanking) or (mapping.qLength - end < flanking):  # segment is clipped given flanking region length
							log.debug('Mapping segment on read {} clipped, not adding'.format(curr_read))
							clipped_segments.append(mapping)
							continue
					new_read_id = curr_read+'.{}'.format(i)
					segment = SeqRecord(new_read_id, reads[curr_read].seq[start-flanking:end+flanking])
					mapping.qName = new_read_id
					new_mappings.append(mapping)

					read_segments.append(segment)
					extracted_reads[curr_read].append(segment)
			else:
				no_valid_segments.append(curr_read)
			curr_read = m.qName
			read_mappings = []
		
		if not is_valid(m):
			continue
		segment = (m.qStart, m.qEnd, m)
		if not read_mappings:
			read_mappings.append(segment)
		if all(x[0] > m.qEnd for x in read_mappings) or all(x[1] < m.qStart for x in read_mappings): # m is up or downstream of all previously seen mapping locations for curr_read
			read_mappings.append(segment)

	log.info('Extracted {} read segments from {} reads'.format(len(read_segments), len(extracted_reads)))
	log.info('{} reads had no valid segments'.format(len(no_valid_segments)))
	log.info('Ignored {} clipped segments'.format(len(clipped_segments)))


	return dict([(x.id, x) for x in read_segments]), new_mappings


if __name__ == '__main__':
	log.level_name = 'WARNING'
	main()