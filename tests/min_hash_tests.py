#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'min_hash_tests.py'

import unittest
import cProfile
from pstats import Stats

from long_read_aligner.kmer import Read
from long_read_aligner.min_hash import MinHash
from long_read_aligner.sketch_directors import GreedyPairSketchDirector, SingleSketchDirector, ExhaustivePairSketchDirector
from long_read_aligner.kmer import KmerSingle, Window
from long_read_aligner.utils.string_utils import CreateWindow, CreateKMERPairs, PerturbWindow


class MinHashTests(unittest.TestCase):

	def test_calcSketch(self):
		"""Test that the Sketch lenght is equal to the number of hashing fucntions"""
		n_hash_permutations = 10
		l_min_hash = MinHash(n_hash_permutations)

		kmers = ['AA', 'TT', 'CC', 'GG']
		kmers = [KmerSingle(i, x) for i, x in enumerate(kmers)]

		Sketch = l_min_hash.calcSketch(kmers)
		self.assertEqual(len(Sketch), n_hash_permutations)

	#def test_getMinimumHash(self):
		# """Test that we can get the minimum hash"""
		# hashes = [0, 1, 2, 3, 4, 5]
		# min_hash = MinHash.getMinimumHash(hashes)
		# self.assertEquals(1, min_hash)

	def test_approximateJaccardGood(self):
		"""Test that the approximate jaccard is 1 at max identity"""
		min_hash = MinHash(3)
		sig1 = [20, 20, 10]
		sig2 = [20, 20, 10]
		jaccard = min_hash.approximateJaccard(sig1, sig2)
		self.assertEquals(1, jaccard)

	def test_approximateJaccardBad(self):
		min_hash = MinHash(3)
		sig1 = [100, 200, 30]
		sig2 = [20, 20, 10]
		jaccard = min_hash.approximateJaccard(sig1, sig2)
		self.assertEquals(0, jaccard)


class SketchDirectorTests(unittest.TestCase):

	def setUp(self):
		"""Enable profiling"""
		self.pr = cProfile.Profile()
		self.pr.enable()
		print "\n<<<---"

	def tearDown(self):
		"""Report profiling results"""
		p = Stats(self.pr)
		p.strip_dirs()
		p.sort_stats('cumtime')
		p.print_stats()
		print "\n--->>>"

	def test_getSketch(self):
		""""""
		test_seq = "AAATTACCCGGG"
		test_seq2 = "AAATTTCCCGGG"
		director = SingleSketchDirector(3, 100, MinHash)
		
		sig1 = director.getSketch(test_seq)
		sig2 = director.getSketch(test_seq2)

		jaccard = director._min_hash.approximateJaccard(sig1, sig2)
		print jaccard
		self.assertTrue(0 < jaccard < 1)

	# def test_time(self):
	# 	"""Kind of a stress test to get some time benchmarks. Here we generate a Sketch for a 1kb window"""
	# 	test_seq = CreateWindow(1000)
	# 	director = SingleSketchDirector(8, 100, MinHash)
	# 	director.getSketch(test_seq)

	def test_getSimilarBig(self):
		"""Kind of a stress test to get some time benchmarks. Here we generate a Sketch for a 1kb window"""
		print "SINGLE"
		test_seq = CreateWindow(1000)
		q_test_seq = PerturbWindow(test_seq, 0.3)
		director = SingleSketchDirector(16, 200, MinHash)
		director.run(test_seq)

		Sketch = director.getSketch(q_test_seq)
		q_test_read = Read(q_test_seq)
		q_test_read.setSketch(Sketch)
		print director.getSimilarReads(q_test_read)

	# def test_SketchCache(self):
	# 	""""""
	# 	print "SINGLE"
	# 	test_seq = "AAATTACCCGGG"
	# 	test_seq2 = "AAATTTCCCGGG"
	# 	test_seq3 = "AAATTTCCAGGT"
	# 	director = SingleSketchDirector(8, 400, MinHash)
		
	# 	director.run(test_seq)
	# 	director.run(test_seq2)
	# 	director.run(test_seq3)
	# 	Sketch = director.getSketch(test_seq)
		
	# 	test_read = Read(test_seq)
	# 	test_read.setSketch(Sketch)
	# 	print director.getSimilarReads(test_read)

	def test_exhaustive(self):
		print "EXHAUSTIVE"
		test_seq = CreateWindow(1000)
		q_test_seq = PerturbWindow(test_seq, 0.3)
		
		director = ExhaustivePairSketchDirector(8, 200, MinHash)
		director.run(test_seq)

		Sketch = director.getSketch(q_test_seq)
		q_test_read = Read(q_test_seq)
		q_test_read.setSketch(Sketch)
		print director.getSimilarReads(q_test_read)

	def test_greedy(self):
		print "GREEDY"
		test_seq = CreateWindow(1000)
		q_test_seq = PerturbWindow(test_seq, 0.3)
		
		director = GreedyPairSketchDirector(8, 200, MinHash)
		director.run(test_seq)

		Sketch = director.getSketch(q_test_seq)
		q_test_read = Read(q_test_seq)
		q_test_read.setSketch(Sketch)
		print director.getSimilarReads(q_test_read)

