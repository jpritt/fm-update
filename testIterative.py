#! /usr/bin/env python
import bwt
import random
import time
import iterativeUpdate as iu 
import sys

''' Compute the number of correctly matched substrands after the original strand is mutated '''

def testIter(genomeLen, numReads, readLen, mutFreq, errors):
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
        base = random.randint(0, len(t2)-1)

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

    return iu.iterativeUpdate(fm, b, alphabet, reads, startsOrig, errors, 5)

# default parameters
genomeLen = 5000    # length of the reference genome
numReads = 500      # number of reads to generate
readLen = 50        # length of generated reads
mutFreq = 0.02      # proportion of bases to mutate
errors = 1          # number of errors to allow in matching

'''
for i in xrange(1,11):
    readLen = 10*i
    print 'Read len: ' + str(readLen)
    numRuns = 5
    avgInit = 0
    avgIter = 0
    for run in xrange(numRuns):
        print '  Run ' + str(run)
        initAcc, iterAcc = testIter(genomeLen, numReads, readLen, mutFreq, errors)
        avgInit += initAcc
        avgIter += iterAcc
    print '  Avg naive accuracy:     ' + str(float(avgInit) / numRuns)
    print '  Avg iterative accuracy: ' + str(float(avgIter) / numRuns)
    sys.stdout.flush()
print '\n'
'''

'''
readLen = 50
for i in xrange(6,11):
    mutFreq = 0.01*i
    print 'Mutation freq: ' + str(mutFreq)
    numRuns = 5
    avgInit = 0
    avgIter = 0
    for run in xrange(numRuns):
        print '  Run ' + str(run)
        initAcc, iterAcc = testIter(genomeLen, numReads, readLen, mutFreq, errors)
        avgInit += initAcc
        avgIter += iterAcc
    print '  Avg naive accuracy:     ' + str(float(avgInit) / numRuns)
    print '  Avg iterative accuracy: ' + str(float(avgIter) / numRuns)
    sys.stdout.flush()
print '\n'
'''

mutFreq = 0.01
numReadsVals = [1000, 1500, 2000]
for numReads in numReadsVals:
    print 'Num Reads: ' + str(numReads)
    numRuns = 5
    avgInit = 0
    avgIter = 0
    for run in xrange(numRuns):
        print '  Run ' + str(run)
        initAcc, iterAcc = testIter(genomeLen, numReads, readLen, mutFreq, errors)
        avgInit += initAcc
        avgIter += iterAcc
    print '  Avg naive accuracy:     ' + str(float(avgInit) / numRuns)
    print '  Avg iterative accuracy: ' + str(float(avgIter) / numRuns)
    sys.stdout.flush()
