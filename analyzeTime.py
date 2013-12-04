#! /usr/bin/env python
import bwt
import random
import time

''' Compare the time to rebuild the index vs updating for substitute, insert, and delete '''

lengths = [100, 1000, 10000, 100000]
runLens = [10, 10, 10, 10]

print 'Size\t\tBuilt\t\t\tInsert\t\t\tDelete\t\t\tSubstitue\t\t\tAverage'

for i in xrange(len(lengths)):
    buildTime = 0
    updateTimes = dict()
    updateTimes[1] = 0
    updateTimes[2] = 0
    updateTimes[0] = 0

    length = lengths[i]
    numRuns = runLens[i]
    for n in xrange(numRuns):
    
        # Generate a long random string of random length
        t = ['$']
        for i in range(length):
            t = [random.choice(['A', 'C', 'T', 'G'])] + t

        alphabet = ['$', 'A', 'C', 'G', 'T']
        b = 50

        # Construct the fm index
        startBuild = time.time()
        fm = bwt.constructFM(t, b, alphabet)
        buildTime += time.time() - startBuild

        letters = set(t)
        letters.remove('$')
    
        # Substitution of a character
        subId = random.randint(0,length-1)
        newChar = random.choice(list(letters))
        t2 = t[:subId] + [newChar] + t[subId+1:]

        startUpdate = time.time()
        fm_new = bwt.substitute(fm, b, alphabet, subId, newChar)
        timeUpdate = time.time() - startUpdate

        updateTimes[0] += timeUpdate


        # Insertion of a character
        insertId = random.randint(0,length-1)
        newChar = random.choice(['A', 'C', 'T', 'G'])
        t2 = t[:insertId] + [newChar] + t[insertId:]
    
        startUpdate = time.time()
        fm_new = bwt.insert(fm, b, alphabet, insertId, newChar)
        timeUpdate = time.time() - startUpdate

        updateTimes[1] += timeUpdate


        # Deletion of a character
        deleteId = random.randint(0,length-1)
        t2 = t[:deleteId] + t[deleteId+1:]
    
        startUpdate = time.time()
        fm_new = bwt.delete(fm, b, alphabet, deleteId)
        timeUpdate = time.time() - startUpdate

        updateTimes[2] += timeUpdate

    buildTime /= numRuns
    for k in updateTimes.keys():
        updateTimes[k] /= numRuns

    print str(length) + '(Update)\t' + str(buildTime) + '\t' + str(updateTimes[1]) + '\t' + str(updateTimes[2]) + '\t' + str(updateTimes[0]) + '\t' + str((updateTimes[1]+updateTimes[2]+updateTimes[0])/3)

