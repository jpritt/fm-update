#! /usr/bin/env python
import bwt
import random
import time
import iterativeUpdateError as iu 
import sys

''' Compute the number of correctly matched substrands after the original strand is mutated '''

def testIter(genomeLen, numReads, readLen, mutFreq, errors, errorFreq):
    # Initialize the reference gene
    t = ['$']
    for i in range(genomeLen):
        t = [random.choice(['A', 'C', 'T', 'G'])] + t

    # constuct the fm index
    alphabet = ['$', 'A', 'C', 'G', 'T']
    b = 5
    fm = bwt.constructFM(t, b, alphabet)

    startsOrig = []
    for i in range(numReads):
        start = random.randint(0, genomeLen-readLen-1)
        startsOrig += [start]
    starts = startsOrig[:]

    # mutate the reference genome to get new genome
    t2 = t[:]
    for i in range(int(round(mutFreq*genomeLen))):
        base = random.randint(2*readLen, len(t2)-2*readLen)

        mutType = random.randint(0,2)
        # substitution
        if mutType == 0:
            t2[base] = random.choice(['A', 'C', 'T', 'G'])
        # insertion
        elif mutType == 1:
            t2 = t2[:base] + [random.choice(['A', 'C', 'T', 'G'])] + t2[base:]
 
            for s in xrange(len(starts)):
                if starts[s] >= base:
                    starts[s] += 1
        # deletion
        else:
            t2 = t2[:base] + t2[base+1:]
            for s in xrange(len(starts)):
                if starts[s] >= base:
                    starts[s] -= 1


    # generate reads from new genome
    reads = []
    for i in xrange(len(starts)):
        reads += [t2[starts[i]:starts[i]+readLen]]
 
        # introduce substitution errors with 1% chance at each base
        for j in xrange(readLen):
            if random.random() < errorFreq:
                reads[i][j] = random.choice(['A', 'C', 'T', 'G'])

    return iu.iterativeUpdateError(fm, b, alphabet, reads, startsOrig, errors, 5, True, readLen, genomeLen)

# default parameters
genomeLen = 5000    # length of the reference genome
numReads = 500      # number of reads to generate
readLen = 50        # length of generated reads
mutFreq = 0.01      # proportion of bases to mutate
errors = 1          # number of errors to allow in matching
errorFreq = 0.01

'''
for i in xrange(1,6):
    numReads = 500*i
    print 'Num reads: ' + str(numReads)
    numRuns = 5
    avgAccuracy = 0
    avgSize = 0
    for run in xrange(numRuns):
        print '  Run ' + str(run)
        accuracy, sizes = testIter(genomeLen, numReads, readLen, mutFreq, errors, errorFreq)
        avgAccuracy += accuracy
        avgSize += sizes[0]
    print '  Avg accuracy: ' + str(float(avgAccuracy) / numRuns)
    print '  Avg size:     ' + str(float(avgSize) / numRuns)
    sys.stdout.flush()
print '\n'
'''

numReads = 1000
for i in xrange(11):
    errorFreq = 0.005*i
    print 'Error freq: ' + str(errorFreq)
    numRuns = 5
    avgAccuracy = 0
    avgSize = 0
    for run in xrange(numRuns):
        print '  Run ' + str(run)
        accuracy, sizes = testIter(genomeLen, numReads, readLen, mutFreq, errors, errorFreq)
        avgAccuracy += accuracy
        avgSize += sizes[0]
    print '  Avg accuracy: ' + str(float(avgAccuracy) / numRuns)
    print '  Avg size:     ' + str(float(avgSize) / numRuns)
    sys.stdout.flush()
print '\n'

