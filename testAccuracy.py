#! /usr/bin/env python
import bwt
import random
import time

''' Compute the number of correctly matched substrands after the original strand is mutated '''

# variable parameters
genomeLen = 5000    # length of the reference genome
numReads = 500      # number of reads to generate
readLen = 50        # length of generated reads
mutFreq = 0.02      # proportion of bases to mutate
errors = 1          # number of errors to allow in matching

def countCorrect(genomeLen, numReads, readLen, mutFreq, errors):
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

    # Match reads against t2
    correct = 0
    incorrect = 0
    for i in range(numReads):
        #print 'Read ' + str(i+1)
        #print '  ' + ''.join(reads[i])
        m = bwt.findApproximate(fm, b, alphabet, ''.join(reads[i]), errors)
        found = False
        #print 'Searching for ' + str(startsOrig[i])
        #print m
        #print ''.join(reads[i])
        #print ''.join(t[starts[i]:starts[i]+readLen])
        for j in xrange(-errors, errors+1):
            if startsOrig[i]+j in m.keys() and not found:
                #print 'Found!\n'
                correct += 1
                found = True
        if not found:
            #print 'Not found\n'
            incorrect += 1
    print '  Accuracy: ' + str(correct) + ' / ' + str(correct+incorrect) + ' = ' + str(float(correct)/(incorrect+correct))
    return float(correct) / (incorrect+correct)


for i in xrange(1,11):
    readLen = 10*i
    print 'Read len: ' + str(readLen)
    numRuns = 5
    avgAccuracy = 0
    for run in xrange(numRuns):
        avgAccuracy += countCorrect(genomeLen, numReads, readLen, mutFreq, errors)
    print '  Avg accuracy: ' + str(float(avgAccuracy) / numRuns)
