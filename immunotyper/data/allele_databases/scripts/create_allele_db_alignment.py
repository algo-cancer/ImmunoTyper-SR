from __future__ import print_function
import sys, os
import argparse
from Bio import SeqIO
from src.common import fasta_from_seq, create_temp_file
from src.tools.msa_wrappers import make_consensus

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
def oprint(*args, **kwargs):
    print(*args, file=sys.stdout, **kwargs)


def convert_imgt_ids(seqs):
    '''Parses IMGT allele sequence ids into a format that clustal allows by removing ( ) / symbols. This prevents id from getting truncated by clustal. '''
    replace_char = lambda x: x.replace('(', '_').replace(')', '_').replace("/", '_')
    
    for x in seqs:
        if 'Novel' not in x.id:
            idx = x.id.split('|')
            allele_id = replace_char(idx[1])
            
            x.id = '|'.join([idx[0]]+[allele_id]+idx[2:])
        else:
            x.id = replace_char(x.id)
            x.seq = x.seq.reverse_complement()
    return seqs

def restore_char(allele_id):
    if 'OR' in allele_id or '_D' in allele_id:
        allele_id = allele_id.replace('_OR', '/OR')
        allele_id = allele_id.replace('_D', '/D')
    if '_I' in allele_id or '_V' in allele_id:
        allele_id = allele_id[:4]+'('+allele_id[5:]
        allele_id = allele_id.replace('_', ')')
    return allele_id


def convert_clustal_ids(seqs):
    '''Parses ids fasta output of Clustal to restore to original IMGT id format'''
    for x in seqs:
        if 'Novel' not in x.id:
            idx = x.id.split('|')
            allele_id = restore_char(idx[1])
            # fix functional 
            functional = idx[3]
            if functional[0] == '_':
                functional = '('+functional[1]+')'
            try:
                x.id = '|'.join([idx[0]]+[allele_id, idx[2], functional]+idx[4:])
            except TypeError as e:
                print((x.id, x.seq))
                raise e
        else:
            x.id = restore_char(x.id)
    return seqs

def main(input_path, output_path=None):
    """Takes the path of a fasta as input, modifys IDs, runs MSA with ClustalW, restores IDs, returns list of SeqIO objects with gapped sequences"""

    # remove whitespaces from fasta ids 
    inp = create_temp_file()
    process = os.system("sed 's/ /_/g' {} > {}".format(input_path, inp.name))

    eprint('Converting ids of sequences')
    imgt = list(SeqIO.parse(inp.name, 'fasta')) 
    converted_imgt = convert_imgt_ids(imgt)
    converted_imgt_file = create_temp_file(write_data=fasta_from_seq(*zip(*[(x.id, x.seq) for x in converted_imgt])),
                                            suffix='.fasta', delete=False)
    
    eprint(f'Running clustal on {len(converted_imgt)} sequences from file {converted_imgt_file.name}')
    output_file = create_temp_file(delete=False)
    if not output_path:
        output_path = output_file.name
    eprint('Running {}'.format('clustalw -INFILE={} -OUTPUT=FASTA -OUTFILE={}'.format(converted_imgt_file.name, output_path)))
    os.system('clustalw -INFILE={} -OUTPUT=FASTA -OUTFILE={}'.format(converted_imgt_file.name, output_path))
    clustal = list(SeqIO.parse(output_path, 'fasta'))
    eprint('Created clustal alignment with {} sequences'.format(len(clustal)))
    converted_imgt_file.close()
    output_file.close()

    return convert_clustal_ids(clustal)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fasta', help="Path to fasta used as input to clustal", required=True)
    args = parser.parse_args()

    output_file = os.path.splitext(os.path.basename(args.input_fasta))[0]+'-aligned.fasta'
    with open(output_file, 'w') as f:
        f.write(fasta_from_seq(*zip(*[(x.id, x.seq) for x in main(args.input_fasta)])))
        eprint('Saving output to {}'.format(output_file))
    alignments = list(SeqIO.parse(output_file, 'fasta'))

    consensus_path = os.path.splitext(os.path.basename(args.input_fasta))[0]+'-consensus.fasta'
    consensus = make_consensus(alignments)

    eprint(f'Saving consensus to {consensus_path}')
    with open(consensus_path, 'w') as f:
        f.write(fasta_from_seq(f'{os.path.splitext(consensus_path)[0]}', consensus))

