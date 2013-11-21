#! /usr/bin/env python
import math
import numpy as np


def rotations(t):
    ''' Return list of rotations of input string t '''
    tt = t * 2
    return [ tt[i:i+len(t)] for i in xrange(0, len(t)) ]

def bwm(t, a):
    ''' Return lexicographically sorted list of t's rotations and suffix array sampled every a indices '''
    rot = rotations(t)
 
    # match each rotation with its index
    rot2 = [(rot[i], i) for i in xrange(len(rot))]
    sort = sorted(rot2)
 
    bwm = [x for (x,_) in sort]
    sa = [y for (_,y) in sort]
  
    return sorted(rotations(t)), sa[(a-1)::a]

def bwtViaBwm(t,a):
    ''' Given T, returns BWT(T) and sampled SA(T) by way of the BWM '''
    bw, sa = bwm(t,a)
    return ''.join(map(lambda x: x[-1], bw)), sa

def rankBwt(bw):
    ''' Given BWT string bw, return parallel list of B-ranks.  Also
        returns tots: map from character to # times it appears. '''
    tots = dict()
    ranks = []
    for c in bw:
        if c not in tots: tots[c] = 0
        ranks.append(tots[c])
        tots[c] += 1
    return ranks, tots

def firstCol(tots):
    ''' Return map from character to the range of rows prefixed by
        the character. '''
    first = {}
    totc = 0
    for c, count in sorted(tots.iteritems()):
        first[c] = (totc, totc + count)
        totc += count
    return first

def reverseBwt(bw):
    ''' Make T from BWT(T) '''
    ranks, tots = rankBwt(bw)
    first = firstCol(tots)
    rowi = 0 # start in first row
    t = '$' # start with rightmost character
    while bw[rowi] != '$':
        c = bw[rowi]
        t = c + t # prepend to answer
        # jump to row that starts with c of same rank
        rowi = first[c][0] + ranks[rowi]
    return t

def getSAIndex(fm, i, a, alphabet):
    (first, last, sa, checkpoints) = fm
    ''' Returns the index in the full SA of the given index i in the bw transform, using the given subsampled SA '''
    if (i+1) % a == 0:
        return sa[(i+1)/a - 1]
    else:
        newi = first[last[i]][0] + getCount(last, i, checkpoints, b, alphabet)
        return (1 + getSAIndex(fm, newi, a, alphabet)) % len(last)

def getCheckpoints(bw, b, alphabet):
    ''' Returns counts of each letter at every b positions in the bw string '''
    counts = dict()
    for l in alphabet:
        counts[l] = 0

    checkpoints = np.zeros((len(bw) / b, len(alphabet)))
    for i in xrange(len(bw)):
        counts[bw[i]] += 1
        if (i % b) == (b-1):
            for k,v in counts.items():
                checkpoints[i/b][alphabet.index(k)] = v
    return checkpoints
   
def getCount(bw, i, checkpoints, b, alphabet):
    '''Find the b-index of the character at the given position in the burrows-wheeler transform''' 
    prevC = i - ((i+1) % b)
    nextC = prevC+b
    if i == prevC:
        return int(checkpoints[(i+1)/b - 1][alphabet.index(bw[i])]) - 1
    elif prevC > -1 and (nextC > len(bw) or (i - prevC) < (nextC - i)):
        count = int(checkpoints[(i+1)/b - 1][alphabet.index(bw[i])])
        for x in bw[prevC+1:i]:
            if x == bw[i]:
                count += 1
        return count
    else:
        count = int(checkpoints[(i+1)/b][alphabet.index(bw[i])])
        for x in bw[(i+1):(nextC+1)]:
            if x == bw[i]:
                count -= 1
        return count-1

def findRange(fm, a, b, substring, start, end):
    (first, last, sa, checkpoints) = fm
    alphabet = sorted(first.keys())

    if len(substring) < 1:
        print 'Error: substring must have length >= 1'
        return []
    else:
        startId = first[substring[-1]][0] + start
        endId = first[substring[-1]][0] + end
        if len(substring) == 1:
            matches = []
            for i in xrange(startId, endId):
                matches.append(getSAIndex(fm, i, a, alphabet))
            return sorted(matches)
        else:
            minId = end-start
            maxId = 0
            for i in xrange(startId, endId):
                if last[i] == substring[-2]:
                    currId  = getCount(last, i, checkpoints, b, alphabet)
                    if currId < minId:
                        minId = currId
                    if currId > maxId:
                        maxId = currId
            if minId > maxId:
                return []
            else:
                return findRange(fm, a, b, substring[:-1], int(minId), int(maxId)+1)

def find(fm, a, b, substring):
    ''' Find all indexes of the substring in the text represented in bw '''
    first = fm[0]
    (minI, maxI) = first[substring[-1]]
    return findRange(fm, a, b, substring, 0, maxI-minI)

def insert(index, c):
    ''' Update the BWT for an insertion of character c into position index in the original string '''


t = 'TTAGCCAATGGAATGGAAGCCGAT$'
query1 = 'AGCC'
query2 = 'AAT'
query3 = 'ACGTA'
alphabet = sorted(set(t))

a = 3
b = 4
last, sa = bwtViaBwm(t, a)
_, tots = rankBwt(last)
first = firstCol(tots)
checkpoints = getCheckpoints(last, b, alphabet)
fm = (first, last, sa, checkpoints)

print 'T = ' + t
print 'Searching for: ' + query1
matches = find(fm, a, b, query1)
print matches

print 'Searching for: ' + query2
matches = find(fm, a, b, query2)
print matches

print 'Searching for: ' + query3
matches = find(fm, a, b, query3)
print matches

