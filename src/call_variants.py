import sys, tempfile, sets
import create_allele_db
from subprocess import Popen, PIPE
import common
from common import log, SeqRecord, fasta_from_seq
from Bio import SeqIO



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
			# print index
			# print seq[:index+1]
			# print forward_map
			# print result
	return (''.join(result), forward_map, backward_map)



def adjust_variant_positions(variants, forward_mapping):
	## For converting variant position from aligned to un-aligned (ie w/o periods) positions
	## Input:	forward_mapping is forward_map output of remove_periods()
	v_changed = []
	for v in variants:
		# new_pos = (forward_mapping[v['pos']] if 'INS' not in v['op'] else forward_mapping[v['pos']]+1)
		new_pos = forward_mapping[v['pos']]
		v_changed.append({'pos': new_pos, 'op': v['op']})
	return v_changed



def get_variants(allele_seq, consensus_seq):
	## Calls variants for allele_seq relative to consensus_seq
	## Input:	Input sequences must be aligned with '.' as the inbdel character
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
	consensus_seq_np, forward, backward = remove_periods(consensus_seq)
	converted_variants = adjust_variant_positions(variants, forward)
	# log.debug(forward)

	msg = msg + test_variants(converted_variants, consensus_seq_np, allele_seq.replace('.', ''), variants)

	return (converted_variants, allele_seq, msg)


