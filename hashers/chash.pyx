#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'hasher.pyx'


cdef class Hasher(object):
	cdef int _a, _b, _c

	def __init__(self, int a, int b):
		self._a = a
		self._b = b
		self._c = max(a, b) - 100

	cpdef long hashIt(self, long x):
		"""
		Hashes the input.
		:param int x: The input value to hash.
		:return int: The hash value.
		"""
		return (x + self._b) % self._c
