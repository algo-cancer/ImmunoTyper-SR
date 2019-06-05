#!/usr/bin/env python
# 786

import os, re, time, sys
import pprint
import logbook, logbook.more

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
	if not (isinstance(name, str)):
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


class SeqRecord(object):
	def __init__(self, id, seq):
		self.id = id
		self.seq = seq

class NamedFunction(object):
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
