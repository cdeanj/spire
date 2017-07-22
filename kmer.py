#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'kmer.py'


import sys
import math
from collections import defaultdict
from abc import ABCMeta, abstractmethod, abstractproperty


class Read(object):

	__slots__ = ['seq', '_sketch']

	def __init__(self, seq):
		self.seq = seq
		self._sketch = None

	def setSketch(self, sketch):
		self._sketch = sketch

	@property
	def sketch(self):
		return self._sketch


class Kmer(object):
	
	__metaclass__ = ABCMeta

	@abstractproperty
	def integer_val(self):
		pass

	@abstractproperty
	def idx(self):
		pass

	@abstractproperty
	def __str__(self):
		pass


class KmerSingle(Kmer):
	
	__slots__ = ['seq', 'idx', 'integer_val']
	
	def __init__(self, idx, seq):
		self.seq = seq
		self.integer_val = hash(seq)
		self.idx = idx

	def __str__(self):
		return self.seq

	def __eq__(self, other):
		return self.seq == other.seq

	def __repr__(self):
		return "<Kmer:{}, idx:{}>".format(self.seq, self.idx)

	def __index__(self):
		return self.idx

	def __gt__(self, other):
		return self.idx > other


class KmerPair(Kmer):

	__slots__ = ['kmer1', 'kmer2', 'distance', 'integer_val']

	def __init__(self, kmer1, kmer2):
		self.kmer1 = kmer1.seq
		self.kmer2 = kmer2.seq
		self.distance = kmer2.idx - kmer1.idx
		self.integer_val = hash(self.kmer1 + self.kmer2 + str(self.distance))

	@property
	def idx(self):
		return 'blah'

	def __str__(self):
		pass

class Window(object):
	
	KMER = KmerSingle

	def __init__(self, seq, kmer_size, n_buckets):
		self.seq = seq
		self.kmer_pairs = None
		self.buckets = getBuckets(n_buckets, len(seq))
		self.kmers = self.makeWords(seq, kmer_size)

	def __getitem__(self, key):
		return self.kmers[key]

	def makeWords(self, w, k):
		return [self.KMER(self.buckets[i], w[i: i + k]) for i in xrange(len(w) - k + 1)]

	def __repr__(self):
		return "<{}:n_kmers:{}>".format(self.__class__.__name__, len(self.kmers))

	def getDownstreamKmers(self, kmer):
		index = self.kmers.index(kmer) + 1
		for kmer in self.kmers[index:]:
			yield kmer

	def makeKmerPairs(self, kmer):
		return [KmerPair(kmer, next_kmer) for next_kmer in self.getDownstreamKmers(kmer)]

def getBuckets(n_buckets, len_seq):
	items_per_bucket = math.ceil(float(len_seq) / float(n_buckets))
	running_bucket = range(n_buckets)
	current_bucket = running_bucket[0]
	buckets = {}
	bucket_count = 0
	
	for i in range(len_seq):
		
		if bucket_count < items_per_bucket:
			
			buckets[i] = running_bucket[current_bucket]
			bucket_count += 1
		
		else:
			current_bucket += 1
			bucket_count = 0
			buckets[i] = running_bucket[current_bucket]
	
	return buckets



