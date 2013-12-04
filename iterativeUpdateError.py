#! /usr/bin/env python
import bwt
import random
import time
import sys
import cPickle

''' Repeatedly match all possible reads to the genome, update, then repeat '''

def iterativeUpdateError(fm, b, alphabet, reads, starts, errors, maxIters, threshold=False, readLen=50, genomeLen=5000, prop=1):
    unmatched = [1]*len(reads)
    numUnmatched = len(reads)
    prevSize = 2*len(reads)
    currIter = 0

    sizes = []
    correct = 0
    incorrect = 0
    while numUnmatched > 0 and currIter < maxIters and float(prevSize - numUnmatched) / prevSize > 0.1:
        #threshold = 0.5 * numUnmatched * readLen / genomeLen
        threshold = 0.25 * numUnmatched * readLen / genomeLen

        currIter += 1
        prevSize = numUnmatched
    
        # Match reads against t2
        mutations = dict()

        # match all reads to genome, collect mutations
        for i in xrange(len(reads)):
            if unmatched[i] == 1:
                m = bwt.findApproximate(fm, b, alphabet, ''.join(reads[i]), errors)
            
                if len(m) > 0:
                    unmatched[i] = 0

                    if i < prop*len(reads):
                        for k,edits in m.items():
                            for v in edits:
                                if v[0] == 2:
                                    vnew = (v[0],v[1]+k)
                                else:
                                    vnew = (v[0],v[1]+k,v[2])
                                if vnew in mutations:
                                    mutations[vnew] += 1
                                else:
                                    mutations[vnew] = 1

                    found = False
                    for j in xrange(-errors, errors+1):
                        if starts[i]+j in m and not found:
                            correct += 1
                            found = True
                    if not found:
                        incorrect += 1

        mutationsString = cPickle.dumps(mutations)
        sizes += [sys.getsizeof(mutationsString)]

        # apply mutations to fm index
        for k,v in mutations.items():
            if v >= threshold:
                if k[0] == 1:
                    fm = bwt.insert(fm, b, alphabet, k[1], k[2])
                    for i in xrange(len(starts)):
                        if starts[i] >= k[1]:
                            starts[i] += 1
                elif k[0] == 0:
                    fm = bwt.substitute(fm, b, alphabet, k[1], k[2])
                elif k[0] == 2:
                    fm = bwt.delete(fm, b, alphabet, k[1])
                    for i in xrange(len(starts)):
                        if starts[i] >= k[1]:
                            starts[i] -= 1
                else:
                    print 'Error: k[0] = ' + str(k[0])
 
        numUnmatched = sum(unmatched)
        print "    Iter " + str(currIter) + " - " + str(correct) + " correct, " + str(incorrect) + " incorrect, " + str(len(reads)-correct-incorrect) + ' unmatched, length = ' + str(len(mutations)) + ', size = ' + str(sys.getsizeof(mutationsString))


    #print "    Accuracy: " + str(float(correct) / len(reads))
    return float(correct) / len(reads), sizes
