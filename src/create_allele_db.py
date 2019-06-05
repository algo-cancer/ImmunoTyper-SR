from Bio import SeqIO
from collections import OrderedDict
import common
from common import log
# from noindent import NoIndent, NoIndentEncoder
import json, os, sets
import call_variants
from subprocess import Popen, PIPE
import Bio




## TODO
# Clean up repetative code in get_variants()
# Remove unnecessary functions (eg fix_truncated_sequences())




def fix_truncated_sequences(allele_seq, consensus_seq):
	## checks if allele sequences is truncated, denoted by '.' at beggining or end of sequnce
	## if yes, replaces truncation '.' by corresponding consensus sequence

	if allele_seq[0] != '.' and allele_seq[-1] != '.':
		return str(allele_seq)
	elif allele_seq[0] == '.':
		for index, pos in enumerate(allele_seq):
			if pos != '.': 
				break
		allele_seq = consensus_seq[:index] + allele_seq[index:]
	else:
		for index, pos in reversed(list(enumerate(allele_seq))):
			if pos != '.':
				break
		allele_seq = allele_seq[:index+1] + consensus_seq[index+1:]
	return str(allele_seq)

def adjust_var_pos(variants, seq):
	## Adjusts variant position to those without the IMGT periods

	nucl_ranges = []
	period_counter = 0

	start = None

	# log.debug(variants)

	## Iterate through seq to find ranges of nucleotides (ie w/o periods) and associated positional shift when periods deleted
	for index, pos in enumerate(seq):


		if start == None and pos != '.':
			# log.debug('Start of range @ pos {}. Start val = {} pos = {}'.format(index, start, pos))
			start = index
		# log.debug('Pos {} Start {} index {} len {}'.format(pos, start, index, len(seq)))
		if start != None and (pos == '.' or index == len(seq)-1):
			nucl_ranges.append((start, (index-1 if index < len(seq)-1 else index), period_counter))
			start = None
		if pos == '.': 
			# log.debug('Incrementing period counter from {} @ pos {}. Context: {} Start = {}'.format(period_counter, index, seq[index-10:index+10], start))
			period_counter += 1
	# log.debug(nucl_ranges)

	for var in variants:
		for start, end, shift in nucl_ranges:
			if var['pos'] >= start and var['pos'] <= end:
				# log.debug('Variant {} is in range {}-{} and will be shifted by {}'.format(var, start, end, shift-1))
				var['pos'] = var['pos'] - shift


def remove_IMGT_periods(seq):
	## Removes IMGT periods from seq value. If not consensus first adjusts variant positions using adjust_var_pos()

	# if 'variants' in list(allele.keys()): adjust_var_pos(allele['variants'], allele['seq'])

	seq = ''.join([x for x in list(seq) if x != '.'])
	return seq

def get_consensus(path=None):
	if not path: path = '../database/V-QUEST-reference-allele-db+no-period-references+consensus.clustalw.consensus-seq.fasta'
	if os.path.exists(path):
		result = SeqIO.parse(path, 'fasta')
		with open(path, 'r') as f: log.debug('Using consensus from {} : \n{}'.format(path, f.read()))
		return str(list(SeqIO.parse(path, 'fasta'))[0].seq)
	else:
		log.error('\nAllele consensus sequence does not exist at {}'.format(path))

def fetch_reference(name):
	
	command = 	['wget -O - -nv',
				'http://www.imgt.org/download/GENE-DB/{}'.format(name),
				"| grep -E 'IGHV.*Homo' -A1 --no-group-separator",
				'> ../database/{}'.format(name)]
	process = Popen(" ".join(command), shell=True, executable='/bin/bash', stdout = PIPE, stderr = PIPE)
	stdout, stderr = process.communicate()

	return True if os.path.exists('../database/{}'.format(name)) else False

def get_conserved_aa(path='../database/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP'):
	if not os.path.exists(path):
		log.info('Fetching IMGT V gene reference amino acid sequences')
		if not fetch_reference('IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP'):
			log.error('Unable to find/fetch from IMGT V gene amino acid references')
			os.sys.exit()

	## iterate through fasta, saving counts of AA at each position in aa_counts
	def add(x, d):
		if x not in d: 
			d[x] = 1
		else: d[x] +=1
	num_references = 0
	aa_counts =[]
	for record in SeqIO.parse(path, 'fasta'):
		num_references += 1
		if not aa_counts:
			for i in range(len(record.seq)):
				aa_counts.append({})
		for index, p in enumerate(str(record.seq)):
			if index > len(aa_counts)-1:
				aa_counts.append({})
			add(p, aa_counts[index])
		

	## get consensus sequence from aa_counts as ('AA', degree of conservation)
	seqAA = []
	for ind, p in enumerate(aa_counts):
		seqAA.append(max(p.iteritems(), key=lambda x: x[1]))

	## get AA w/ > conservation in > 94% of alleles
	## position stored as nucleotide NOT AA position
	cons_aa = [(i*3, x[0]) for i, x in enumerate(seqAA) if x[1] >= int(0.94*num_references) and x[0] != '.'] 

	## get codon conversion table
	rev_translation = {}
	for k, v in Bio.Data.CodonTable.standard_dna_table.forward_table.iteritems():
		if v in rev_translation:
			rev_translation[v].append(k.lower())
		else:
			rev_translation[v] = [k.lower()]

	## convert positions to position indicies for consensus sequence with periods removed
	## since that will be used for read alignment
	consnp, forward, backward = call_variants.remove_periods(get_consensus())

	## get conserved codons and return
	return dict([(forward[x[0]], rev_translation[x[1]]) for x in cons_aa])



def iterate_fasta(filename='../database/V-QUEST-reference-allele-db+no-period-references.clustalw.fasta'):
	def strip_fasta_ID(s): # strips out allele name from the IMGT fasta naming convention
		s = s.split('|')
		return (s[1], s[3], s[0]) # (name, functional_value, accession)

	if not os.path.exists(filename):
		log.info('Fetching IMGT V gene reference nucleotide sequences')
		if not fetch_reference('IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'):
			log.error('Unable to find/fetch from IMGT V gene nucleotide references')
			os.sys.exit()

	consensus = ''
	current_gene = None

	return_db = OrderedDict()
	consensus = get_consensus().lower().replace('-', '.')

	for record in SeqIO.parse(filename, 'fasta'):
		allele, functional, accession = strip_fasta_ID(record.description)

		# if allele != 'IGHV1-18*01':
		# 	continue 

		# log.debug('\nAllele {}'.format(allele))

		seq = str(record.seq).lower().replace('-', '.')
		# log.debug('Allele: '+ allele)



		if not current_gene or allele.split('*')[0] != current_gene: # first entry in iteration
			current_gene = allele.split('*')[0]
			# return_db[current_gene] = OrderedDict({'alleles': OrderedDict()})
		# log.debug('allele seq:\n{}\nConsensus seq:\n{}'.format(seq, consensus))
		length = len(seq.replace('.', ''))
		variants, seq, msg = call_variants.get_variants(seq, consensus)
		# log.debug('Has {} variants'.format(len(variants)))
		if msg:
			log.warn('\nAllele {}'.format(allele))
			common.log_msg(msg)

				
		return_db[allele] = OrderedDict({	'imgt_accession': accession,
																	'functional': functional,
																	'seq': seq,
																	'length': length})
		if variants:
			# log.debug(variants)
			if __name__ != '__main__':
				variants = sets.Set([(x['pos'], x['op']) for x in variants])
			return_db[allele]['variants'] = variants
	# remove_IMGT_periods(return_db[current_gene]['alleles'][allele])


	return return_db

def prettify_json_output(allele_db):
	for alleles in allele_db.values():
		if 'variants' in list(alleles.keys()):
			alleles['variants'] =  [NoIndent(x) for x in alleles['variants']]
	log.info('Saving database')
	return allele_db


# if __name__ == '__main__':
# 	log.level_name = 'WARNING'
# 	with open('output_get_variants_from_ref_fasta_test.yml', 'w') as f:
# 			print >>f, json.dumps(prettify_json_output(iterate_fasta('../database/V-QUEST-reference-allele-db.fasta')), indent=2, cls=NoIndentEncoder)
