import sys
import os
import read
import assemble23


if __name__ == '__main__':

    with open(sys.argv[1], "r") as f:   
        K, startContigsNumber, sim, maxDeep, startLen, minMatch, scoreThreshold, supportNumberThreshold, coverRatio = read.read_parameter(f)
    
    assemble23.find_best_path(sys.argv[2], sim, startLen, scoreThreshold, supportNumberThreshold)
