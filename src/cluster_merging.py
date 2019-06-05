from common import log, SeqRecord
from collections import OrderedDict
from minimap_wrapper import MinimapWrapper

class ReMapper(object):

	remapped_clusters = {}

	def __init__(self, to_remap, mapper=MinimapWrapper(params="-cx map-pb -k10 -w3 -N5")):
		self.mapper=mapper
		if to_remap:
			self.to_remap = to_remap
		else:
			self.to_remap=lambda x: True if len(x) ==1 else False


	def merge_clusters(self, clusters):
		log.info('Remapping clusters')

		self.remapped_clusters = {}

		querys = {}
		targets = OrderedDict()
		for c in clusters:
			if self.to_remap(c):
				log.debug('Remapping reads from cluster {}'.format(c.id))
				for r in c:
					querys[r.id] = r
					r.original_cluster=c
			else:
				targets[str(c.id)] = c

		log.info('Remapping {} reads'.format(len(querys)))

		query_seqs = querys.values()
		target_seqs = [SeqRecord(c.id, c.consensus.replace('.', '')) for c in targets.values()]
		
		mapping = self.mapper.run(query_seqs, target_seqs)


		## sorted list insertion function
		def add_mapping(new, arr, less_than=lambda x, y: True if x.total_errors() < y.total_errors() else False):
			if less_than(new, arr[0]):
				return [new] + arr
			if len(arr) == 1 or not less_than(new, arr[-1]):
				return arr + [new]
			for i, m in enumerate(arr):
				if less_than(new, m):
					return arr[:i] + [new] + arr[i:]

		sort = {}
		for m in mapping:
			m.quality = m.total_errors()
			try:
				sort[m.qName] = add_mapping(m, sort[m.qName])
			except KeyError:
				sort[m.qName] = [m]

		for r_id, maps in sort.iteritems():
			try:
				read = querys[r_id]
				read.cluster_mappings = maps
				targets[maps[0].tName].reads[read.id] = read 		## add to mapped cluster
				read.cluster = targets[maps[0].tName]
				targets[maps[0].tName].has_remaps = True
				if maps[0].tName not in self.remapped_clusters:
					self.remapped_clusters[maps[0].tName] = [read]
				else:
					self.remapped_clusters[maps[0].tName].append(read)
			except KeyError as e:
				print(maps[0])
				print(maps[0].tName in targets)
				print(sorted(targets.keys()))
				print(targets[maps[0].tName].reads.keys())
				raise e
		return targets.values()



class RemapToParent(ReMapper):

	def merge_clusters(self, clusters):

		for c in clusters:
			c.has_remaps = False
		clusters = super(RemapToParent, self).merge_clusters(clusters)

		remapped_clusters = [c for c in clusters if c.has_remaps]

		result_clusters = []
		to_be_removed = []
		for c in remapped_clusters:
			parent = c.parent
			if parent in result_clusters: # this parent has already been a merge target, do not add again
				continue
			all_sibling_reads = c.reads
			to_be_removed.append(c)
			for i in clusters:
				if i.parent == parent and i.id != c.id: # i is a sibling of c
					to_be_removed.append(i)
					all_sibling_reads.update(i.reads)

			parent.reads = all_sibling_reads
			result_clusters.append(parent)

		result_clusters.extend([c for c in clusters if c not in to_be_removed])

		return result_clusters




