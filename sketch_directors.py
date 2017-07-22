__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = "min_hash.py"


import sys
import numpy as np
import operator
from collections import defaultdict, Counter
from abc import ABCMeta, abstractmethod

from hasher import HashFactory
from kmer import Kmer, Window, Read

import pyximport; pyximport.install()
from hashers.chash import Hasher



class SketchDirector(object):
	"""Flow control for generating minhash Sketchs for a given sequence"""

	__metaclass__ = ABCMeta

	def __init__(self, word_size, n_hashing_functions, hasher):
		"""
		Init method for class
		:param int word_size: The kmer size to use
		:param n_hashing_functions: The number of hashing functions to use.
		"""
		self.word_size = word_size
		self._n_hashing_functions = n_hashing_functions
		self._min_hash = hasher(n_hashing_functions)
		self._Sketch_cache = {x:defaultdict(list) for x in xrange(n_hashing_functions)}

	@abstractmethod
	def getSketch(self):
		pass

	def run(self, seq):
		Sketch = self.getSketch(seq)
		self.storeSketch(seq, Sketch)
		return Sketch

	def storeSketch(self, seq, Sketch):
		for i, min_mer in enumerate(Sketch):
			self._Sketch_cache[i][min_mer].append(seq)

	def querryRead(self, Sketch):
		for i, min_mer in enumerate(Sketch):
			yield self._Sketch_cache[i][min_mer]

	def getSimilarReads(self, read):
		
		similar_reads = []
		for q_read_list in self.querryRead(read.sketch):
			for q_read in q_read_list:
				if q_read != read.seq:
					similar_reads.append(q_read)
		
		return Counter(similar_reads)


class GreedyPairSketchDirector(SketchDirector):


	def __init__(self, word_size, n_hashing_functions, hasher):
		super(GreedyPairSketchDirector, self).__init__(word_size, n_hashing_functions, hasher)
	
	def getSketch(self, seq):
		"""
		Gets the minhash Sketch for a given sequence
		:param str seq: A DNA sequence
		:return list(int): A list representing the minhash Sketch
		"""
		
		window = Window(seq, self.word_size, 2)
		
		primary_Sketch = self._min_hash.calcSketch(window.kmers)
		kmer_pairs = []
		
		for min_hash in primary_Sketch:
			kmer_pairs.extend(window.makeKmerPairs(min_hash[1]))
		
		return [x[0] for x in self._min_hash.calcSketch(kmer_pairs)]

	@classmethod
	def canCalculate(cls, opts):
		return opts.mode == "greedy"


class SingleSketchDirector(SketchDirector):
	
	def __init__(self, word_size, n_hashing_functions, hasher):
		super(SingleSketchDirector, self).__init__(word_size, n_hashing_functions, hasher)

	def getSketch(self, seq):
		"""
		Gets the minhash Sketch for a given sequence
		:param str seq: A DNA sequence
		:return list(int): A list representing the minhash Sketch
		"""
		window = Window(seq, self.word_size, 2)
		return [x[0] for x in self._min_hash.calcSketch(window.kmers)]

	@classmethod
	def canCalculate(cls, opts):
		return opts.mode == "fast"


class ExhaustivePairSketchDirector(SketchDirector):

	def __init__(self, word_size, n_hashing_functions, hasher):
		super(ExhaustivePairSketchDirector, self).__init__(word_size, n_hashing_functions, hasher)


	def getSketch(self, seq):
		"""
		Gets the minhash Sketch for a given sequence
		:param str seq: A DNA sequence
		:return list(int): A list representing the minhash Sketch
		"""
		
		window = Window(seq, self.word_size, 2)
		
	
		pairs = [pair for kmer in window.kmers for pair in window.makeKmerPairs(kmer)]
		return [x[0] for x in self._min_hash.calcSketch(pairs)]

	@classmethod
	def canCalculate(cls, opts):
		return opts.mode == "deep"