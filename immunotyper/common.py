#!/usr/bin/env python
# 786

import os, sys
import pkg_resources
import logbook, logbook.more
import tempfile
from Bio.Seq import Seq

DATA_DIR = 'immunotyper.data'
def resource_path(key, data_dir_path=DATA_DIR):
    """
    Args:
        key (str): resource to be extracted from the ``data'' directory.
        data_dir_path (str): Relative path from the package root directory.
    Returns:
        str: Full path of the resource.
    """
    return pkg_resources.resource_filename(data_dir_path, key)

db_resource_path = lambda x: resource_path(x, 'immunotyper.data.allele_databases')
databases = {'ighv': {'db_fasta_path': db_resource_path('IMGT_IGHV_reference_allele_db-updated+oscar_novel-aligned+orphon+Ns-aligned.fasta'),
                                                                'consensus_path': db_resource_path('IMGT_IGHV_reference_allele_db-updated+oscar_novel-aligned+orphon+Ns-consensus.fasta'),
                                                                'ignored_alleles_path': resource_path('ignored_alleles.txt')},
                'iglv': {'db_fasta_path': db_resource_path('IGLV-IMGT-allele-db-aligned.fasta'),
                                                                'consensus_path': db_resource_path('IGLV-IMGT-allele-db-consensus.fasta')},
                'trav': {'db_fasta_path': db_resource_path('IMGT_TRAV_reference_allele_db-aligned.fasta'),
                                                            'consensus_path': db_resource_path('IMGT_TRAV_reference_allele_db-consensus.fasta')},
                'trbv': {'db_fasta_path': db_resource_path('TRBV-IMGT-allele-db-aligned.fasta'),
                                                            'consensus_path': db_resource_path('TRBV-IMGT-allele-db-consensus.fasta')},
                'trdv': {'db_fasta_path': db_resource_path('TRDV-IMGT-allele-db-aligned.fasta'),
                                                            'consensus_path': db_resource_path('TRDV-IMGT-allele-db-consensus.fasta')},
                'trgv': {'db_fasta_path': db_resource_path('TRGV-IMGT-allele-db-aligned.fasta'),
                                                            'consensus_path': db_resource_path('TRGV-IMGT-allele-db-consensus.fasta')},
                'igkv':  {'db_fasta_path': db_resource_path('IGKV-IMGT-allele-db-aligned.fasta'),
                                                            'consensus_path': db_resource_path('IGKV-IMGT-allele-db-consensus.fasta')}}



allele_db_mapping_path = {'ighv': db_resource_path('IMGT_IGHV_reference_allele_db-updated+oscar_novel+orphon+Ns.fa'),
                        'iglv': db_resource_path('IGLV-IMGT-allele-db-no_duplicates+Ns.fasta'),
                        'trav': db_resource_path('IMGT_TRAV_reference_allele_db+Ns.fasta'),
                        'trbv': db_resource_path('TRBV-IMGT-allele-db-no_duplicates+Ns.fa'),
                        'trdv': db_resource_path('TRDV-IMGT-allele-db+Ns.fa'),
                        'trgv': db_resource_path('TRGV-IMGT-allele-db+Ns.fa'),
                        'igkv': db_resource_path('IGKV-IMGT-allele-db-no_duplicates+Ns.fa')}


def header(string):
    return '\n\n' + '-'*len(string) + string + '-'*len(string) + '\n\n'

def colorize(text, color='green'):
    return logbook._termcolors.colorize(color, text)

log = logbook.Logger('Cypiripi')
# log.level = logbook.DEBUG

def initialize_logger(debug_log_path='immunotyper-debug'):
    LOG_FORMAT = '{record.message}'
    if debug_log_path:
        debug_log_path = debug_log_path+'.log'
        if os.path.exists(debug_log_path):
            os.remove(debug_log_path)
        
        handler = logbook.NestedSetup([logbook.NullHandler(),
                                    logbook.FileHandler(debug_log_path, level='DEBUG', format_string=LOG_FORMAT),
                                    logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level='INFO', bubble=True)])
    else:
        handler = logbook.NestedSetup([logbook.NullHandler(),
                                    logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level='INFO', bubble=True)])

    handler.push_application()



def log_msg(msg):
    for m in msg:
        f = m[0]
        f(m[1])

mean = lambda x: round(float(sum(x))/len(x), 3)
median = lambda x: sorted(x)[len(x)/2]
def mean_median(x):
    return '{} mean, {} median'.format(mean(x), median(x))


def fasta_from_seq(name, seq):
    ## Input: Name, seq can be str or iterable yielding str
    result = []
    if not isinstance(name, str):
        try:
            for n, s in zip(name, seq):
                result.append('>{}\n{}'.format(n, s))
        except TypeError:
            log.error('Please provide a iterable or string')
            raise TypeError
    else:
        result = ['>{}\n{}'.format(name, seq)]
    return '\n'.join(result)


def get_columns(data):
    col_width = max(len(str(word)) for row in data for word in row) + 2  # padding
    result = []
    for row in data:
        result.append("".join(str(word).ljust(col_width) for word in row))
    return '\n'.join(result)

def print_columns(data):
    print(get_columns(data))


class Read():
    n_buffer = 100 # Number of flanking Ns added to allele database sequences
    def __init__(self, iden, seq):
        self.id = iden
        self.seq = seq
        self.mappings = []
        self.mappings_dict = {}
        self.primary_mapping = None

    def __getitem__(self, key):
        if isinstance(key, slice):
            indices = range(*key.indices(len(self.seq)))
            return ''.join(self.seq[i] for i in indices)
        return self.seq[key]
    
    def __str__(self):
        return str(self.seq)
    
    def __len__(self):
        return len(self.seq)
    
    def __hash__(self):
        return hash(self.id)
    
    def __eq__(self, other):
        return self.id == other.id

    @property
    def reverse_complement(self):
        return str(Seq(self.seq).reverse_complement())
    
    @property
    def is_discarded(self):
        return True if self.discard_var.X > 0 else False
    
    def add_mapping(self, mapping):
        self.mappings.append(mapping)
        self.mappings_dict[self.reference_id_parser(mapping.reference_name)] = mapping
        try:
            self.primary_mapping = self.get_primary_mapping(self.mappings)
        except ValueError:
            pass

    def get_distance(self, allele):
        '''Gets NM tag edit distance from mapping to provided allele'''
        if isinstance(allele, str):
            allele_id = allele
        else:
            allele_id = allele.id
        
        m = self.mappings_dict[allele_id]
        num_Ns = abs(min(m.reference_start-self.n_buffer, 0)) + abs(min(len(self.allele_db[allele_id])-m.reference_end+self.n_buffer, 0))

        return m.get_tag('NM')-num_Ns

    def get_start(self, allele):
        '''Gets reference_start from mapping to provided allele'''

        if isinstance(allele, str):
            allele_id = allele
        else:
            allele_id = allele.id
        return max(0, self.mappings_dict[allele_id].reference_start-self.n_buffer)
    
    def get_unclipped_start(self, allele):
        '''Gets reference_start from mapping to provided allele'''
        #TODO test
        start = self.get_start(allele)
        end = self.get_end(allele)
        length_diff = len(self.seq) - (end-start)
        if length_diff == 0:
            return start
        else:
            return start-length_diff

    def get_end(self, allele):
        '''Gets reference_end from mapping to provided allele'''
        if isinstance(allele, str):
            allele_id = allele
        else:
            allele_id = allele.id
        return min(len(self.allele_db[allele_id]), self.mappings_dict[allele_id].reference_end-self.n_buffer)
    
    def get_unclipped_end(self, allele):
        '''Gets reference_start from mapping to provided allele'''
        #TODO test
        start = self.get_start(allele)
        end = self.get_end(allele)
        length_diff = len(self.seq) - (end-start)
        if length_diff == 0:
            return end
        else:
            return end+length_diff
    
    def covers_position(self, pos, allele):
        '''Returns true if mapping to allele covers pos'''
        return True if (self.get_start(allele) <= pos and pos <= self.get_end(allele)) else False

    @staticmethod
    def reference_id_parser(raw_id):
        if '|' in raw_id:
            return raw_id.split('|')[1]
        else:
            return raw_id

    @staticmethod
    def get_primary_mapping(mapping_list):
        '''Takes a list of pysam.AlignedSegment objects, returns the one that is primary mapping
        Raises ValueError if it does not have,or has multiple, primary mappings'''
        primary_mapping = [x for x in mapping_list if not x.is_secondary and not x.is_supplementary]
        if len(primary_mapping) > 1:
            s = "\n".join([str(x)+f'\t Reference = {x.reference_name}' for x in primary_mapping])
            log.error(f'More than 1 primary mapping:\n{s}')
        if len(mapping_list) and not primary_mapping:
            raise ValueError('Provided mappings have no primary mapping'.format(str([str(x) for x in mapping_list])))
        return None if not primary_mapping else primary_mapping[0]


class ReadFlanking(Read):
    """For Reads with mapping to allele database that contains flanking sequences. Overides many of the mapping functions to accomodation flanking sequences"""
    def __init__(self, allele_db, *args, **kwargs):
        self.allele_db = allele_db
        self.coding_distances = {}
        super().__init__(*args, **kwargs)
    
    def get_distance(self, allele, coding_only=True):
        '''Gets NM tag edit distance from mapping to provided allele'''
        if isinstance(allele, str):
            allele_id = allele
        else:
            allele_id = allele.id
        if coding_only and allele_id in self.coding_distances:
            return self.coding_distances[allele_id]
        else:
            return self.mappings_dict[allele_id].get_tag('NM')

    def get_start(self, allele, coding_only:bool=False):
        '''Gets reference_start from mapping to provided allele
        coding_only (bool) if True returns relative to coding start (negative = before start)'''

        if isinstance(allele, str):
            allele_id = allele
        else:
            allele_id = allele.id

        if coding_only: 
            try:
                return self.mappings_dict[allele_id].reference_start - self.allele_db[allele_id].coding_start
            except AttributeError:
                pass
        return self.mappings_dict[allele_id].reference_start
    
    def get_end(self, allele, coding_only:bool=False):
        '''Gets reference_end from mapping to provided allele
                coding_only (bool) if True returns relative to coding end (negative = after end)'''

        if isinstance(allele, str):
            allele_id = allele
        else:
            allele_id = allele.id
        if coding_only: 
            try:
                return self.allele_db[allele_id].coding_end - self.mappings_dict[allele_id].reference_end
            except AttributeError:
                pass
        return self.mappings_dict[allele_id].reference_end

    def covers_position(self, pos, allele, coding_only:bool=True):
        '''Returns true if mapping to allele covers pos
            coding_only (bool) if True position is relative to start of coding aka doesnt factor in flanking'''
        if self.allele_db[allele].has_flanking and len(self.allele_db[allele].upstream_flanking) > 0:
            pos = len(self.allele_db[allele].upstream_flanking) + pos
        return True if (self.get_start(allele, coding_only=False) <= pos and pos <= self.get_end(allele, coding_only=False)) else False



class SeqRecord(Read):
    '''Implemented for backwards compatibility'''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    

class NamedFunction():
    def __init__(self, name, f):
        self.f = f
        self.name = name

    def __call__(self, *args, **kwargs):
        return self.f(*args, **kwargs)

    def __str__(self):
        return self.name


class PrintRedirect:
    def __init__(self, output_path=None):
        if not output_path: 
            self.output_path = os.devnull
        else:
            self.output_path = output_path
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(self.output_path, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def create_temp_file(write_data=None, suffix='', delete=True, mode="w", **kwargs):
    '''Use fasta_from_seq(*zip(*[(x.id, x.seq) for x in SeqRecords]) to write seqrecord sequences
    Use tempfile.name to get path'''
    result = tempfile.NamedTemporaryFile(suffix=suffix, delete=delete, mode=mode, **kwargs)
    if write_data:
        result.seek(0)
        result.write(write_data)
        result.truncate()
        result.flush()
    return result

from contextlib import contextmanager

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

# import subprocess
# import shlex
# def run_command(command_string: str, check: bool=True, log_stout: bool=False, log_stderr: bool=False, dont_split: bool=False, *args, **kwargs):
#         '''Wrapper for subprocess.run
#         Args
#             command_string
#             check                       run check arg value
#             log_stout, log_sterr        writes stout, sterr to log.ino if True respectively
#             *args, **kwargs             Additional args for run
#         Raises subprocess.CalledProcessError, writes stdout, sterr to log.error if check=True and non-zero
#             exit status returned '''
#         try:
#             command_string = shlex.split(command_string) if not dont_split else command_string
#             result = subprocess.run(command_string, check=check, capture_output=True, *args, **kwargs)
#         except subprocess.CalledProcessError as e:
#             log.error(f"Command {command_string} returned non-zero exit status:\n{e.stdout.decode('UTF-8')}\n{e.stderr.decode('UTF-8')}")
#             raise e
        
#         if log_stout: log.info(result.stdout)
#         if log_stderr: log.info(result.stderr)

def run_command(command_string: str, *args, **kwargs):
    os.system(command_string)
