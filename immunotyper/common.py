#!/usr/bin/env python
# 786

import os, sys
import pkg_resources
import logbook, logbook.more
import tempfile
from Bio.Seq import Seq
import subprocess

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
# ... existing code ...

class AlleleDatabaseMissingException(Exception):
    """Exception raised when the allele database for a given gene type is missing."""
    pass

supported_gene_types = {'IGHV', 'IGKV', 'IGLV', 'TRBV', 'TRDV', 'TRGV', 'TRAV'}
def get_database_config(gene_type):
    """
    Returns the database configuration for the specified gene type.

    Args:
        gene_type (str): The gene type for which to retrieve the database configuration.

    Returns:
        dict: A dictionary containing paths to the database files for the specified gene type.

    Raises:
        AlleleDatabaseMissingException: If the gene type is not supported.
    """
    if gene_type.upper() not in supported_gene_types:
        raise AlleleDatabaseMissingException(f"Allele database for gene type '{gene_type}' is missing.")

    base_config = {
        'db_fasta_path': '{}/{}-IMGT-allele-db-aligned.fasta',
        'consensus_path': '{}/{}-IMGT-allele-db-consensus.fasta',
        'gene_clusters_path': '{}/{}-IMGT-allele-db-gene_clusters.tsv',
        'gene_type': gene_type.upper()
    }

    # Update the paths with the specific gene type
    for key, template in base_config.items():
        if template and '{}' in template:
            base_config[key] = db_resource_path(template.format(gene_type.upper(), gene_type.upper()))

    # Add ignored_alleles_path for IGHV specifically
    if gene_type.lower() == 'ighv':
        base_config['ignored_alleles_path'] = db_resource_path('IGHV/IGHV-ignored_alleles.txt')

    return base_config

def get_allele_db_mapping_path(gene_type):
    """
    Returns the database resource path for the specified gene type.

    Args:
        gene_type (str): The gene type for which to retrieve the database path.

    Returns:
        str: The full path of the database resource.

    Raises:
        AlleleDatabaseMissingException: If the gene type is not supported.
    """
    if gene_type.upper() not in supported_gene_types:
        raise AlleleDatabaseMissingException(f"Allele database for gene type '{gene_type}' is missing.")

    return db_resource_path(f'{gene_type.upper()}/{gene_type.upper()}-IMGT-allele-db-no_duplicates+Ns.fa')

# ... existing code ...

def header(string):
    return '\n\n' + '-'*len(string) + string + '-'*len(string) + '\n\n'

def colorize(text, color='green'):
    return logbook._termcolors.colorize(color, text)

# Create the logger but don't set handlers yet
log = logbook.Logger('ImmunoTyper')

# Add a basic NullHandler to start
basic_handler = logbook.NullHandler()
basic_handler.push_application()

def initialize_logger(debug_log_path='immunotyper-debug'):
    """Initialize the logger with appropriate handlers.
    
    Args:
        debug_log_path: Path for debug log file. If None, only console logging is enabled.
    """
    LOG_FORMAT = '{record.message}'
    
    # Create new handler setup
    handlers = [logbook.NullHandler()]
    
    if debug_log_path:
        debug_log_path = debug_log_path + '.log'
        if os.path.exists(debug_log_path):
            os.remove(debug_log_path)
            
        handlers.extend([
            logbook.FileHandler(debug_log_path, level='DEBUG', format_string=LOG_FORMAT),
            logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level='INFO', bubble=True)
        ])
    else:
        handlers.append(
            logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level='INFO', bubble=True)
        )

    # Setup and push the new handler configuration
    handler_setup = logbook.NestedSetup(handlers)
    handler_setup.push_application()



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

def run_command(command: str, check: bool = True, quiet: bool = False) -> None:
    """
    Run a shell command with optional output suppression.
    
    Args:
        command: The shell command to execute
        check: If True, raises RuntimeError on command failure
        quiet: If True, suppresses stdout and stderr during normal execution
    
    Raises:
        RuntimeError: If the command fails and check is True
    """
    log.debug(f'Running {command}')
    try:
        # Configure stdout/stderr handling
        stdout = subprocess.DEVNULL if quiet else None
        stderr = subprocess.PIPE
        
        # Run the command using subprocess.run()
        result = subprocess.run(command, 
                              shell=True, 
                              check=check, 
                              stdout=stdout,
                              stderr=stderr,
                              text=True)
    except subprocess.CalledProcessError as e:
        # If the command fails and check is True, raise an exception with stderr
        stderr_output = e.stderr.strip()
        raise RuntimeError(f"Command '{command}' failed with error: {str(e)}\nStderr: {stderr_output}")
