

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = "to the moon"

import argparse
import sys
from min_hash import MinHash
from sketch_directors import GreedyPairSketchDirector, SingleSketchDirector, ExhaustivePairSketchDirector
from utils.string_utils import PerturbWindow

def parse_cmdline_params(cmdline_params):
    
    info = "Calculate a Jaccard index for a querry sequence against a reference"
    
    parser = argparse.ArgumentParser(description=info)
    
    parser.add_argument('-q', '--querry_string', type=str, required=True,
                        help='Provide a querry string')

    parser.add_argument('-r', '--reference_sequence', type=str, required=True,
    					help='Please provide a reference string')

    parser.add_argument('-p', '--perturb_input', type=float, 
    					help='probabilty of a sequencing error at any given base')

    parser.add_argument('-n', '--n_hashing_functions', type=int, default=100, 
    					help='Number of hashing f(x) to considering when calculating a signature')

    parser.add_argument('-k', '--word_size', type=int, help="kmer_size")

    parser.add_argument('-m', '--mode', type=str, choices = ['greedy', 'fast', 'deep'], 
    					default='fast', help="The algorithm to generate the kmer-pair scheme")

    return parser.parse_args(cmdline_params)

if __name__ == '__main__':

	modes = [GreedyPairSketchDirector, SingleSketchDirector, ExhaustivePairSketchDirector]
	opts= parse_cmdline_params(sys.argv[1:])
	
	for mode in modes:
		if mode.canCalculate(opts):
			director = mode(opts.word_size, opts.n_hashing_functions, MinHash)
			break
		else:
			continue

	sketch1 = director.run(opts.reference_sequence)
	sketch2 = director.run(opts.querry_string)
	jaccard = director._min_hash.approximateJaccard(sketch1, sketch2)
	print jaccard
