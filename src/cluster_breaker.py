from clustering import EditDistanceCalculator
from common import log, fasta_from_seq, SeqRecord, get_columns, PrintRedirect
from cluster_class import SubCluster
from collections import Counter
import abc, os, sys
import msa_generators, lpinterface
from msa_generators import Poa
import call_variants
from minimap_wrapper import MinimapWrapper
from lpinterface import BreakerModel
from Bio.Seq import Seq

class ClusterBreaker(object):
	__metaclass__ = abc.ABCMeta

	def __init__(self, to_break=None):
		if to_break:
			self.to_break = to_break
		else:
			self.to_break = lambda x: x 	# filter function to determine which clusters to break



	@abc.abstractmethod
	def break_cluster(self, cluster):
		raise NotImplementedError

	def break_clusters(self, clusters, *args, **kwargs):
		log.debug('Breaking clusters with:\n{}'.format(str(self)))
		result = [] 
		for i, cluster in enumerate(clusters):
			if self.to_break(cluster):
				try:
					sub_clusters = self.break_cluster(cluster, *args, **kwargs)
					if not sub_clusters:
						log.warn('Cluster {} not broken'.format(cluster.id))
						result.append(cluster)
					else:
						log.info('Breaking cluster {} into {} sub_clusters'.format(cluster.id, len(sub_clusters)))
						result.extend(sub_clusters)
				except (lpinterface.NoSolutionsError, UnboundLocalError, TypeError, ValueError) as e:
					log.error('Cluster breaking failed for cluster {} - see log'.format(cluster.id))
					log.debug(sys.exc_info())
					result.append(cluster)
			else:
				result.append(cluster)
		return result




class IlpBreaker(ClusterBreaker):
	"""
	Args
		poa_cache (func):			Takes cluster as argument, returns string path to save / load location of poa_cache
	"""
	top_candidate_filter = 10

	def __init__(self, allele_db,
						read_depth,
						max_copy=4,
						minimap=MinimapWrapper(params='-cx map-pb -k10 -w3'),
						poa='../tools/poaV2/poa', 
						poa_cache=None,
						to_break=None,
						gurobi_log_path=None):
		self.allele_db = allele_db
		self.read_depth = read_depth
		self.max_copy=max_copy
		self.min_var_support = int(read_depth*0.9)
		self.max_var_support = int(read_depth*2)
		self.min_candidate_support = 0.5*read_depth

		self.minimap=minimap
		self.poa = msa_generators.Poa(src=poa)
		self.poa_cache = poa_cache
		self.gurobi_log_path=gurobi_log_path

		log.debug(str(self))

		super(IlpBreaker, self).__init__(to_break)

	def __str__(self):
		result = 'IlpBreaker\n'
		parameters = [('Read depth:', self.read_depth),
						('Min variant support:', self.min_var_support),
						('Max variant support:', self.max_var_support),
						('Top candidate filter:', self.top_candidate_filter),
						('Min candidate support:', self.min_candidate_support)]
		result = result+get_columns(parameters)
		return result

	def set_model_parameters(self, cluster, model):
		""" Sets default model parameters"""

		model.maxnum=int(round(len(cluster)/(self.read_depth*0.9)))
		model.minnum=int(round(len(cluster)/(self.read_depth*1.1)))
		model.minsize=int(round(len(cluster)/(self.read_depth*0.9)))
		model.expcov=int(self.read_depth)
		model.maxcopy = self.max_copy

	def __str__(self):
		result = 'ILP Cluster breaker'
		return result

	def break_cluster(self, cluster, model=None, *args, **kwargs):
		"""
		Args
			*args, **kwargs: 			Overrides default BreakerModel parameters (set by self.set_model_parameters)
		"""
		try:
			if not cluster.consensus:
				raise AttributeError
		except AttributeError:
			raise ValueError('Cluster {} does not have consensus'.format(cluster.id))

		log.info('-----\nBreaking cluster {}\n-----'.format(cluster.id))
		log.info('Size: {}'.format(len(cluster)))
		if self.poa_cache:
			try:
				cache_dir = os.path.dirname(self.poa_cache(cluster))
				if not os.path.exists(cache_dir):
					log.debug('Cache path does not exist, creating {}'.format(cache_dir))
					os.makedirs(cache_dir)
				self.poa.parse_poa_pir(self.poa_cache(cluster), cluster)
				log.info('Using poa cache: {}'.format(self.poa_cache(cluster)))
			except (KeyError, IOError):
				self.poa.generate_consensus_msa(cluster, output_path=self.poa_cache(cluster))	
		else:
			self.poa.generate_consensus_msa(cluster)


		self.get_read_variants(cluster)

		self.get_candidates(cluster)

		self.get_candidate_variants(cluster)

		self.filter_variants(cluster)

		for r in cluster:
			r.candidate = None

		if not model:
			model = BreakerModel()
			self.set_model_parameters(cluster, model)
		model.build(cluster, *args, **kwargs)

		log.debug('Running model...')
		with PrintRedirect(self.gurobi_log_path):
			model.solve()
		
		return self.get_subclusters(cluster, model)
		
	def get_subclusters(self, cluster, model):

		calls, results = model.get_assignments()
		sub_clusters = []
		for i, (call, reads) in enumerate([(call, reads) for call, reads in calls if reads]):
			sub = SubCluster([cluster.reads[r_id] for r_id in reads],
							cluster_id=i,
							parent=cluster,
							clustering_tool=str(self))
			sub.breaker_call = [call]
			sub.breaker_results = results
			sub.model = model
			sub_clusters.append(sub)
		return sub_clusters
	


	def warn_candidates(self, cluster):
		from common import get_columns
		ground_truth = [x.id.split('_')[0] for x in cluster]

		count = get_columns(Counter(ground_truth).most_common())

		missing = set(ground_truth).difference(set([x.id for x in cluster.candidates]))
		if missing:
			log.warn('---\nCluster {}\nMissing {} from candidates\nGround Truth:\n{}'.format(cluster.id, str(missing), count))

	def get_candidates(self, cluster):
		from collections import Counter
		from common import SeqRecord
		log.debug('Getting candidates...')

		# # get new mappings for the reads using the filtered allele database
		# reads = {}
		# for r in cluster:
		# 	r.mapping = []
		# 	reads[r.id] = r
		# mappings = self.minimap.run(reads.values(),
		# 				'../database/Complete.Human.IGHV_IMGT.Feb2018.Corey.filtered.fasta',
		# 				params='-cx map-pb -k10 -w3 -N5')
		# for m in mappings:
		# 	reads[m.qName].mapping.append(m)
		# missed_from_mapping = set(reads.keys()).difference(set([x.qName for x in mappings]))
		# if missed_from_mapping:
		# 	log.warn('{} reads did not map to any alleles'.format(len(missed_from_mapping)))
		# 	log.debug('\n'.join(missed_from_mapping))

		## let candidates be all ambiguous top mappings, for every read in cluster
		candidates = []
		for r in cluster:
			maps = sorted(r.mapping, key=lambda x: x.tErrors(), reverse=True)
			r.candidates=[]
			for m in maps:
				if m.tErrors() <= maps[0].tErrors(): 
					candidates.append(m.tName.split('|')[1])
					r.candidates.append(m.tName.split('|')[1])
		
		## Now filter based on count
		candidates = Counter(candidates).most_common()
		log.debug('{} candidate before filter'.format(len(candidates)))
		# log.debug(get_columns(candidates))

		candidates_filtered = set([x for x in candidates[:min([len(candidates), self.top_candidate_filter])] if x[1] > self.min_candidate_support])
		log.info('{} candidates after top {} filter'.format(len(candidates_filtered), self.top_candidate_filter))

		if len(candidates_filtered) == 0:
			raise ValueError('Cluster {} has no candidates after filtering\n{}'.format(cluster.id, str(candidates[:min([len(candidates), self.top_candidate_filter])])))

		cluster.candidates = [SeqRecord(x[0], self.allele_db[x[0]]['seq'].replace('.', '')) for x in candidates_filtered]
		for cand, support in zip(cluster.candidates, candidates_filtered):
			cand.read_mapping_support = support[1]

		self.fix_candidate_rev_comp(cluster)




	def fix_candidate_rev_comp(self, cluster):
		mappings = self.minimap.run([SeqRecord('cons', cluster.consensus.replace('.', ''))],
									cluster.candidates,
									params=self.minimap.params + ' -N{}'.format(len(cluster.candidates)))

		mappings = dict([(m.tName, m) for m in mappings])

		for c in cluster.candidates:
			try:
				if mappings[c.id].strand == '-':
					c.seq = str(Seq(c.seq).reverse_complement())
					log.debug('Reverse complemented candidate {}'.format(c.id))
			except KeyError:
				log.warn('Breaking cluster {}: candidate {} not in rev-compl mapping. Read support: {}'.format(cluster.id, c.id, c.read_mapping_support))


	def get_candidate_variants(self, cluster):
		input_seqs = cluster.candidates + [SeqRecord('cons', cluster.consensus.replace('.', ''))]
		msa = self.poa.generate_msa(input_seqs, 
									params = '-read_fasta {} -preserve_seqorder -pir {} ../tools/poaV2/blosum80.mat')
		cluster.candidate_msa=msa
		cons = msa['cons']
		for c in cluster.candidates:
			var = call_variants.get_variants(msa[c.id], cons)[0]
			c.variants = set(['{}.{}'.format(x['pos'], x['op']) for x in var if 'SNP' in x['op']])

	def filter_variants(self, cluster):
		## filter candidate variants
		all_cand_var = set()
		not_useful = cluster.candidates[0].variants
		for c in cluster.candidates:
			all_cand_var = all_cand_var.union(c.variants)
			not_useful = not_useful.intersection(c.variants)
		log.debug('{} candidate variants before filtering'.format(len(all_cand_var)))
		all_cand_var = set()
		for c in cluster.candidates:
			c.raw_variants = c.variants
			c.variants = c.variants.difference(not_useful)
			all_cand_var = all_cand_var.union(c.variants)
		log.debug('{} candidate variants after filtering - {} removed'.format(len(all_cand_var), len(not_useful)))

		## filter read variants
		all_read_var = []
		for r in cluster:
			all_read_var.extend(list(r.variants))
		log.debug('{} read variants before filtering'.format(len(all_read_var)))
		all_read_var = Counter(all_read_var).most_common()

		## filter read variants based on read support
		valid_var = set([x[0] for x in all_read_var if self.max_var_support >= x[1] and x[1] >= self.min_var_support])
		## ensure all candidate variants stay in there
		valid_var = valid_var.union(all_cand_var)
		log.debug('{} total valid read variants after filtering'.format(len(valid_var)-len(all_cand_var)))

		for r in cluster:
			r.variants = r.variants.intersection(valid_var)



	def get_read_variants(self, cluster):
		for r in cluster:
			var = call_variants.get_variants(r.alignment, cluster.consensus)[0]
			r.variants = set(['{}.{}'.format(x['pos'], x['op']) for x in var if 'SNP' in x['op']])