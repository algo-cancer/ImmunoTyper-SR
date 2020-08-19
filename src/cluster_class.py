from collections import Sequence, OrderedDict, Counter
from common import get_columns, print_columns

class Cluster(Sequence):	
	reads = None			## Orderdict {read_id, SeqIO.SeqRecord-like object}
	aligned_reads = None	##
	clustering_tool = None
	consensus = None		## str of consensus sequence
	consensus_method = None	## str description of consensus method - if using consensus_builder class just set as str(consensus_builder_instance)
	call = None				## Allele call for cluster based on candidate filtering criteria
	candidates = None		## list of candidates
	candidates_method = None ## str description of candidate calling method - if using allele_caller class just set as str(consensus_builder_instance)
	cluster_mappings= []

	def __init__(self, reads, cluster_id, clustering_tool=None):

		self.reads = OrderedDict([(x.id, x) for x in reads])
		if len(reads) != len(self.reads):
			raise ValueError("Cluster reads must have unique id's")
		self.id = cluster_id

		self.seqs = (x.seq for x in self)

		if clustering_tool:
			self.clustering_tool = clustering_tool

		super(Cluster, self).__init__()


	def __getitem__(self, i):
		return list(self.reads.values())[i]
	def __len__(self):
		return len(self.reads)

	def __str__(self):
		result = []
		result.append('Cluster reads:')
		result.extend([x.id for x in self])

		if self.clustering_tool:
			result.append(str(self.clustering_tool))

		if self.consensus:
			result.append('Has consensus{}{}'.format(
														(' and MSA of reads ' if self.aligned_reads else ''), 
														(' created with\n'+ self.consensus_method if self.consensus_method else '')
													)
						)
		if self.call:
			result.append('Allele call: '+ str(self.call))
		if self.candidates:
			result.append('Top {} candidates:\n{}'.format(min(3, len(self.candidates)),
															'\n'.join([str(x) for x in self.candidates[:min(3, len(self.candidates))]])))
		# result.append(self.print_remaps())

		return '\n'.join(result)

	def set_call(self, call):
		if not call:
			self.call = None
		else:
			if isinstance(call, str):
				self.call = [Call(call, self)]
			else:
				self.call = [Call(x, self) for x in call]


	def print_call_summary(self):
		result = []
		if self.call:
			result.append('TOP CALL: \n'+ str(self.call))
		else:
			result.append('NO CALL - '+('no consensus' if not self.consensus else 'no candidate mappings'))

		if self.candidates:
			result.append('\nTop {} candidates:\n'.format(min(3, len(self.candidates))))
			try:
				result.append(str(self.candidates[0].get_header()))
			except AttributeError as e:
				pass
			result.append('\n'.join([str(x) for x in self.candidates[:min(3, len(self.candidates))]]))

		return '\n'.join(result)

	def print_remaps(self, all_attr=False):
		result = []

		data = self.get_remaps(all_attr=all_attr)
		for r, maps in data.iteritems():
			result.append('---\n'+str(r.id))
			for m in maps:
				result.append('{}\t{}'.format(m[0], m[1]) if not all_attr else str(m))

		if result:
			if all_attr:
				result = ['Cluster id\tQuality'] + result
			result = ['\nCluster remappings'] + result
		return '\n'.join(result)

	def get_remaps(self, all_attr=False):
		result = {} 
		for r_id, r in self.reads.iteritems():
			if 'cluster_mappings' not in dir(r):
				continue
			result[r_id] = []
			for m in r.cluster_mappings:
				try:
					if not all_attr:
						result[r_id].append((m.tName, m.quality))
					else:
						result[r_id].append(m)
				except AttributeError as e:
					log.error('Read {} has cluster mappings but quality values have not been set'.format(r.id))
					log.debug(str(m))
					raise e
		return result





	def copy(self):
		new = Cluster(list(self), self.id)

		for a in ['aligned_reads','clustering_tool','consensus','consensus_method','call','candidates', 'candidates_method', 'cluster_mappings']:
			exec('new.{}=self.{}'.format(a, a))
		return new



class SubCluster(Cluster):
	"""docstring for SubCluster"""
	def __init__(self, reads, cluster_id, parent, clustering_tool=None):
		self.parent = parent
		cluster_id = '{}.{}'.format(parent.id, cluster_id)
		super(SubCluster, self).__init__(reads, cluster_id, clustering_tool)



class Call(str):
## str wrapper that has .cluster attr that points to associate cluster
	def __new__(cls, call, cluster):
		obj = str.__new__(cls, call)
		obj.cluster = cluster
		return obj

