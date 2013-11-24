#! /usr/bin/env python
import math
import numpy as np
import random

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
#    sa_complete = [y for (_,y) in sort]
    sa = dict()
    for i in xrange(0, len(sort), a):
        sa[i] = sort[i][1]
  
    return sorted(rotations(t)), sa

def bwtViaBwm(t,a):
    ''' Given T, returns BWT(T) and sampled SA(T) by way of the BWM '''
    bw, sa = bwm(t,a)
    return map(lambda x: x[-1], bw), sa

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
    if i in sa:
        return sa[i]
    else:
        newi = first[last[i]][0] + getCount(last, i, checkpoints, b, alphabet)
        return (1 + getSAIndex(fm, newi, a, alphabet)) % len(last)

def getRowBySA(fm, index, a, alphabet):
    ''' Return the row  in the bwt with the given SA value '''
    (first, last, sa, checkpoints) = fm

    minK = 0
    minDiff = len(last)
    for k,v in sa.items():
        if v == index:
            return k
        elif (v-index) > 0 and (v-index) < minDiff:
            minK = k
            minDiff = v-index

    row = minK
    for i in range(minDiff):
       row = LF(fm, row)
    return row

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
    elif prevC > -1 and (nextC >= len(bw) or (i - prevC) < (nextC - i)):
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
            minId = len(last)
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

def moveRow(fm, i, j):
    ''' Moves the row from i to j, shifting all rows in-between and updates SA and checkpoints '''
    (first, last, sa, checkpoints) = fm
    if i > j:
        last = last[:j] + [last[i]] + last[j:i] + last[i+1:]

        # Update SA
        for k,v in sa.items():
            if k >= j and k < i:
                sa[k+1] = v
                del sa[k]
            elif k == i:
                sa[j] = v
                del sa[k]

        # Loop through checkpoints between j and i
        indexAdded = alphabet.index(last[j])
        for x in xrange(j + b - (j%b+1), i, b):
            checkpoints[(x+1) / b - 1][indexAdded] += 1
            checkpoints[(x+1) / b - 1][alphabet.index(last[x+1])] -= 1
    else:
        last = last[:i] + last[i+1:j+1] + [last[i]] + last[j+1:]
        
        # Update SA
        for k,v in sa.items():
            if k > i and k <= j:
                sa[k-1] = v
                del sa[k]
            elif k == i:
                sa[j] = v
                del sa[k]

        # Loop through checkpoints between j and i
        indexRemoved = alphabet.index(last[j])
        for x in xrange(i + b - (i%b+1), j, b):
            checkpoints[(x+1) / b - 1][indexRemoved] -= 1
            checkpoints[(x+1) / b - 1][alphabet.index(last[x])] += 1
    return (first, last, sa, checkpoints)

def find(fm, a, b, substring):
    ''' Find all indexes of the substring in the text represented in bw '''
    first = fm[0]
    (minI, maxI) = first[substring[-1]]
    return findRange(fm, a, b, substring, 0, maxI-minI)

def insert(fm, index, c):
    ''' Update the BWT for an insertion of character c into position index in the original string '''
    (first, last, sa, checkpoints) = fm

    # Update old row
    row = getRowBySA(fm, index, a, alphabet)
    tempC = last[row]
    last[row] = c

    # Update first column
    for k in alphabet:
        if k == c:
            first[k] = (first[k][0], first[k][1]+1)
        elif k > c:
            first[k] = (first[k][0]+1, first[k][1]+1)

    # Add new row
    bVal = getCount(last, row, checkpoints, b, alphabet)
    if not tempC == c:
        bVal += 1

    newRow = first[c][0] + bVal
    last = last[:newRow] + [tempC] + last[newRow:]


    # Update SA
    for k,v in sa.items():
        if k >= newRow:
            if v >= index:
                v += 1
            del sa[k]
            sa[k+1] = v
        else:
            if v >= index:
                sa[k] = v+1

    print checkpoints
    # Add new checkpoint row
    if len(last) % b == 0:
        # WHY ISN'T THIS CLONING
        newCheckpoint = checkpoints[len(checkpoints)-1][:]
        for x in xrange(len(checkpoints)*b,len(last)):
            newCheckpoint[alphabet.index(last[x])] += 1
        checkpoints = np.vstack([checkpoints, newCheckpoint])
        print checkpoints
    # Update checkpoints
    if row < newRow:
        # Loop through checkpoints between row and newRow
        indexAdded = alphabet.index(c)
        indexRemoved = alphabet.index(tempC)
        for x in xrange(row + b - (row%b+1), newRow, b):
            checkpoints[(x+1) / b - 1][indexAdded] += 1
            checkpoints[(x+1) / b - 1][indexRemoved] -= 1
        for x in xrange(newRow + b - (newRow%b+1), len(last), b):
            checkpoints[(x+1) / b - 1][indexAdded] += 1
            checkpoints[(x+1) / b - 1][alphabet.index(last[x+1])] -= 1
    else:
        # Loop through checkpoints between row and newRow
        indexAdded = alphabet.index(c)
        indexRemoved = alphabet.index(tempC)
        for x in xrange(newRow + b - (newRow%b+1), row, b):
            checkpoints[(x+1) / b - 1][indexRemoved] += 1
            checkpoints[(x+1) / b - 1][alphabet.index(last[x+1])] -= 1
        for x in xrange((row+1) + b - ((row+1)%b+1), len(last), b):
            checkpoints[(x+1) / b - 1][indexAdded] += 1
            checkpoints[(x+1) / b - 1][alphabet.index(last[x+1])] -= 1
    print checkpoints

    

    # Rearrange rows that are now out of order
    fm = (first, last, sa, checkpoints)

    for i in xrange(len(last)):
        if ((i+1) % b == 0):
            n = (i+1)/b-1
            print last[i] + '\t' + str(getSAIndex(fm, i, a, alphabet)) + '\t' + str(checkpoints[n][0]) + '\t' + str(checkpoints[n][1]) + '\t' + str(checkpoints[n][2]) + '\t' + str(checkpoints[n][3]) + '\t' + str(checkpoints[n][4])
        else:
            print last[i] + '\t' + str(getSAIndex(fm, i, a, alphabet))
    print sa
    print ''

    j = getRowBySA(fm, index-1, a, alphabet)
    j2 = LF(fm, newRow)
    while not j == j2:
        newJ = LF(fm, j)
        fm = moveRow(fm, j, j2)
        j = newJ
        j2 = LF(fm, j2)

    return fm 

def LF(fm, index):
    ''' Step forward one step in the bwt and return the next row '''
    (first, last, sa, checkpoints) = fm

    c = last[index]
    return (first[c][0] + getCount(last, index, checkpoints, b, sorted(first.keys()))) % len(last)

length = 10
t = ['$']
for i in range(length):
    t = [random.choice(['A', 'C', 'T', 'G'])] + t
print t

alphabet = sorted(set(t))
a = 3
b = 3
last, sa = bwtViaBwm(t, a)
_, tots = rankBwt(last)
first = firstCol(tots)
checkpoints = getCheckpoints(last, b, alphabet)
fm = (first, last, sa, checkpoints)


for i in xrange(len(last)):
    print last[i] + '\t' + str(getSAIndex(fm, i, a, alphabet))
print '\n'

insertId = random.randint(0,20)
newChar = random.choice(['A', 'C', 'T', 'G'])
print str(insertId) + ', ' + newChar
t2 = t[:insertId] + [newChar] + t[insertId:]
print t2

last2, sa2 = bwtViaBwm(t2, a)
_, tots2 = rankBwt(last2)
first2 = firstCol(tots2)
checkpoints2 = getCheckpoints(last2, b, alphabet)
fm2 = (first2, last2, sa2, checkpoints2)

fm_new = insert(fm, insertId, newChar)

for i in xrange(len(last2)):
    print last2[i] + '\t' + str(getSAIndex(fm2, i, a, alphabet)) + '\t' + fm_new[1][i] + '\t' + str(getSAIndex(fm_new, i, a, alphabet))



'''
t = 'TTAGCCAATGGAATGGAAGCCGAT$'
query1 = 'AGCC'
query2 = 'AAT'
query3 = 'ACTGG'
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



t2 = 'TTAGCCAACTGGAATGGAAGCCGAT$'

last2, sa2 = bwtViaBwm(t2, a)
_, tots2 = rankBwt(last2)
first2 = firstCol(tots2)
checkpoints2 = getCheckpoints(last2, b, alphabet)
fm2 = (first2, last2, sa2, checkpoints2)

print 'T = ' + t2
print 'Searching for: ' + query1
matches = find(fm2, a, b, query1)
print matches

print 'Searching for: ' + query2
matches = find(fm2, a, b, query2)
print matches

print 'Searching for: ' + query3
matches = find(fm2, a, b, query3)
print matches


print '\nUpdating initial index'
fm_updated = insert(fm, 8, 'C')
print 'T = ' + t2
print 'Searching for: ' + query1
matches = find(fm_updated, a, b, query1)
print matches

print 'Searching for: ' + query2
matches = find(fm_updated, a, b, query2)
print matches

print 'Searching for: ' + query3
matches = find(fm_updated, a, b, query3)
print matches


for i in xrange(len(last)):
    for c in alphabet:
        if i >= first[c][0] and i < first[c][1]:
            print c + '\t' + last[i] + '\t' + str(getSAIndex(fm,i,a,alphabet)) + '\t' + str(getCount(last, i, checkpoints, b, alphabet))

print '\n'
(first, last, _, checkpoints) = fm2
for i in xrange(len(last)):
    for c in alphabet:
        if i >= first[c][0] and i < first[c][1]:
            print c + '\t' + last[i] + '\t' + str(getSAIndex(fm2, i, a, alphabet)) + '\t' + str(getCount(last, i, checkpoints, b, alphabet))
'''
