#! /usr/bin/env python
import bwt
import random
import time
import iterativeEM
import iterativeEMDist
import iterativeUpdateError as iu 
import sys
import copy

''' Compute the number of correctly matched substrands after the original strand is mutated '''

def iterEM(genomeLen, numReads, readLen, mutFreq, errors, errorFreq):
    ''' Prop = proportion of reads to contribute to mutations '''
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
 
        # introduce substitution errors with 1% chance at each base
        for j in xrange(readLen):
            if random.random() < errorFreq:
                reads[i][j] = random.choice(['A', 'C', 'T', 'G'])

    tempReads = reads[:]
    tempFM = copy.deepcopy(fm)
    tempStarts = startsOrig[:]
    accuracyOrig, sizeOrig = iterativeEM.iterativeEM(tempFM, b, alphabet, tempReads, tempStarts, errors, 5, readLen, genomeLen, 1)

    tempFM = copy.deepcopy(fm)
    tempReads = reads[:]
    tempStarts = startsOrig[:]
    accuracyRed, sizeRed = iterativeEM.iterativeEM(tempFM, b, alphabet, tempReads, tempStarts, errors, 5, readLen, genomeLen, 0.1)
    
    tempFM = copy.deepcopy(fm)
    tempReads = reads[:]
    tempStarts = startsOrig[:]
    accuracyDist1, sizeDist1 = iterativeEMDist.iterativeEMDist(tempFM, b, alphabet, tempReads, tempStarts, errors, 5, readLen, genomeLen, 1, 1)

    tempFM = copy.deepcopy(fm)
    tempReads = reads[:]
    tempStarts = startsOrig[:]
    accuracyDist2, sizeDist2 = iterativeEMDist.iterativeEMDist(tempFM, b, alphabet, tempReads, tempStarts, errors, 5, readLen, genomeLen, 2, 1)

    tempFM = copy.deepcopy(fm)
    tempReads = reads[:]
    tempStarts = startsOrig[:]
    accuracyDistChunk2, sizeDistChunk2 = iterativeEMDist.iterativeEMDist(tempFM, b, alphabet, tempReads, tempStarts, errors, 5, readLen, genomeLen, 2, 50)
    return (accuracyOrig, accuracyRed, accuracyDist1, accuracyDist2, accuracyDistChunk2), (sizeOrig, sizeRed, sizeDist1, sizeDist2, sizeDistChunk2)

# default parameters
genomeLen = 5000    # length of the reference genome
numReads = 500      # number of reads to generate
readLen = 50        # length of generated reads
mutFreq = 0.01      # proportion of bases to mutate
errors = 1          # number of errors to allow in matching
errorFreq = 0.01

labels = '\tOrig\tReduced\tDist Depth=1\tDist Depth=2\tChunk Depth=2'

for i in xrange(1,9):
    readLen = 10*i
    print 'Read Len: ' + str(readLen)
    numRuns = 5
    avgAccuracies = [0,0,0,0,0]
    avgSizes = [0,0,0,0,0] 
    for run in xrange(numRuns):
        print '  Run ' + str(run)
        accuracies, sizes = iterEM(genomeLen, numReads, readLen, mutFreq, errors, errorFreq)
        for i in xrange(5):
            avgAccuracies[i] += float(accuracies[i]) / numRuns
            avgSizes[i] += float(sizes[i]) / numRuns
        sys.stdout.flush()
    
    print labels
    print '  Acc:\t' + str(avgAccuracies[0]) + '\t' + str(avgAccuracies[1]) + '\t' + str(avgAccuracies[2]) + '\t\t' + str(avgAccuracies[3]) + '\t\t' + str(avgAccuracies[4])
    print '  Size:\t' + str(avgSizes[0]) + '\t' + str(avgSizes[1]) + '\t' + str(avgSizes[2]) + '\t\t' + str(avgSizes[3]) + '\t\t' + str(avgSizes[4])
    sys.stdout.flush()
print '\n'

readLen = 50
for i in xrange(1,6):
    numReads = 500*i
    print 'Num Reads: ' + str(numReads)
    numRuns = 5
    avgAccuracies = [0,0,0,0,0]
    avgSizes = [0,0,0,0,0] 
    for run in xrange(numRuns):
        print '  Run ' + str(run)
        accuracies, sizes = iterEM(genomeLen, numReads, readLen, mutFreq, errors, errorFreq)
        for i in xrange(5):
            avgAccuracies[i] += float(accuracies[i]) / numRuns
            avgSizes[i] += float(sizes[i]) / numRuns
        sys.stdout.flush()
    
    print labels
    print '  Acc:\t' + str(avgAccuracies[0]) + '\t' + str(avgAccuracies[1]) + '\t' + str(avgAccuracies[2]) + '\t\t' + str(avgAccuracies[3]) + '\t\t' + str(avgAccuracies[4])
    print '  Size:\t' + str(avgSizes[0]) + '\t' + str(avgSizes[1]) + '\t' + str(avgSizes[2]) + '\t\t' + str(avgSizes[3]) + '\t\t' + str(avgSizes[4])
    sys.stdout.flush()
print '\n'

numReads = 1000
for i in xrange(11):
    errorFreq = 0.005*i
    print 'Error freq: ' + str(errorFreq)
    numRuns = 5
    avgAccuracies = [0,0,0,0,0]
    avgSizes = [0,0,0,0,0] 
    for run in xrange(numRuns):
        print '  Run ' + str(run)
        accuracies, sizes = iterEM(genomeLen, numReads, readLen, mutFreq, errors, errorFreq)
        for i in xrange(5):
            avgAccuracies[i] += float(accuracies[i]) / numRuns
            avgSizes[i] += float(sizes[i]) / numRuns
        sys.stdout.flush()
    
    print labels
    print '  Acc:\t' + str(avgAccuracies[0]) + '\t' + str(avgAccuracies[1]) + '\t' + str(avgAccuracies[2]) + '\t\t' + str(avgAccuracies[3]) + '\t\t' + str(avgAccuracies[4])
    print '  Size:\t' + str(avgSizes[0]) + '\t' + str(avgSizes[1]) + '\t' + str(avgSizes[2]) + '\t\t' + str(avgSizes[3]) + '\t\t' + str(avgSizes[4])
    sys.stdout.flush()
print '\n'

