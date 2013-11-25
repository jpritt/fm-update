#! /usr/bin/env python
import math
import numpy as np
import random

def rotations(t):
    ''' Return list of rotations of input string t '''
    tt = t * 2
    return [ tt[i:i+len(t)] for i in xrange(0, len(t)) ]

def bwm(t):
    ''' Return lexicographically sorted list of t's rotations and suffix array sampled every a indices '''
    rot = rotations(t)
 
    # match each rotation with its index
    rot2 = [(rot[i], i) for i in xrange(len(rot))]
    sort = sorted(rot2)
 
    bwm = [x for (x,_) in sort]
    sa = [y for (_,y) in sort]
    #sa = dict()
    #for i in xrange(0, len(sort), a):
    #    sa[i] = sort[i][1]
  
    return sorted(rotations(t)), sa

def bwtViaBwm(t):
    ''' Given T, returns BWT(T) and sampled SA(T) by way of the BWM '''
    bw, sa = bwm(t)
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


#def getSAIndex(fm, i, a, alphabet):
#    (first, last, sa, checkpoints) = fm
#    ''' Returns the index in the full SA of the given index i in the bw transform, using the given subsampled SA '''
#    if i in sa:
#        return sa[i]
#    else:
#        newi = first[last[i]][0] + getCount(last, i, checkpoints, b, alphabet)
#        return (1 + getSAIndex(fm, newi, a, alphabet)) % len(last)
#
#def getRowBySA(fm, index, alphabet):
#    ''' Return the row  in the bwt with the given SA value '''
#    (first, last, sa, checkpoints) = fm
#
#    minK = 0
#    minDiff = len(last)
#    for k,v in sa.items():
#        if v == index:
#            return k
#        elif (v-index) > 0 and (v-index) < minDiff:
#            minK = k
#            minDiff = v-index
#
#    row = minK
#    for i in range(minDiff):
#       row = LF(fm, row)
#    return row


def getCheckpoints(bw, b):
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
   
def getCount(bw, i, checkpoints, b):
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
        for x in bw[i:(nextC+1)]:
            if x == bw[i]:
                count -= 1
        return count

def findRange(fm, b, substring, start, end):
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
                matches.append(sa[i])
            return sorted(matches)
        else:
            minId = len(last)
            maxId = 0
            for i in xrange(startId, endId):
                if last[i] == substring[-2]:
                    currId  = getCount(last, i, checkpoints, b)
                    if currId < minId:
                        minId = currId
                    if currId > maxId:
                        maxId = currId
            if minId > maxId:
                return []
            else:
                return findRange(fm, b, substring[:-1], int(minId), int(maxId)+1)

def moveRow(fm, i, j):
    ''' Moves the row from i to j, shifting all rows in-between and updates SA and checkpoints '''
    (first, last, sa, checkpoints) = fm
    if i > j:
        last = last[:j] + [last[i]] + last[j:i] + last[i+1:]

        # Update SA
        sa = sa[:j] + [sa[i]] + sa[j:i] + sa[i+1:]

        # Loop through checkpoints between j and i
        indexAdded = alphabet.index(last[j])
        for x in xrange(j + b - (j%b+1), i, b):
            checkpoints[(x+1) / b - 1][indexAdded] += 1
            checkpoints[(x+1) / b - 1][alphabet.index(last[x+1])] -= 1
    else:
        last = last[:i] + last[i+1:j+1] + [last[i]] + last[j+1:]
        
        # Update SA
        sa = sa[:i] + sa[i+1:j+1] + [sa[i]] + sa[j+1:]

        # Loop through checkpoints between j and i
        indexRemoved = alphabet.index(last[j])
        for x in xrange(i + b - (i%b+1), j, b):
            checkpoints[(x+1) / b - 1][indexRemoved] -= 1
            checkpoints[(x+1) / b - 1][alphabet.index(last[x])] += 1
    return (first, last, sa, checkpoints)

def find(fm, b, substring):
    ''' Find all indexes of the substring in the text represented in bw '''
    first = fm[0]
    (minI, maxI) = first[substring[-1]]
    return findRange(fm, b, substring, 0, maxI-minI)

def insert(fm, index, c):
    ''' Update the BWT for an insertion of character c into position index in the original string '''
    (first, last, sa, checkpoints) = fm

    # Update old row
    row = sa.index(index)
    tempC = last[row]
    last[row] = c

    # Update first column
    for k in alphabet:
        if k == c:
            first[k] = (first[k][0], first[k][1]+1)
        elif k > c:
            first[k] = (first[k][0]+1, first[k][1]+1)

    # update checkpoints        
    indexAdded = alphabet.index(c)
    indexRemoved = alphabet.index(tempC)
    for x in xrange(row + b - (row%b+1), len(last), b):
        checkpoints[(x+1) / b - 1][indexAdded] += 1
        checkpoints[(x+1) / b - 1][indexRemoved] -= 1


    # Add new row
    bVal = getCount(last, row, checkpoints, b)

    newRow = first[c][0] + bVal
    last = last[:newRow] + [tempC] + last[newRow:]


    # Update SA
    for i in xrange(len(sa)):
        if sa[i] >= index:
            sa[i] += 1
    sa = sa[:newRow] + [index] + sa[newRow:]

    # Update checkpoints
    indexRemoved = alphabet.index(tempC)
    for x in xrange(newRow + b - (newRow%b+1), len(last)-1, b):
        checkpoints[(x+1) / b - 1][indexRemoved] += 1
        checkpoints[(x+1) / b - 1][alphabet.index(last[x+1])] -= 1
    # Add new checkpoint row
    if len(last) % b == 0:
        newCheckpoint = np.copy(checkpoints[-1])
        for x in xrange(len(checkpoints)*b,len(last)):
            newCheckpoint[alphabet.index(last[x])] += 1
        #newCheckpoint[alphabet.index(tempC)] += 1
        checkpoints = np.vstack([checkpoints, newCheckpoint])

    # Rearrange rows that are now out of order
    fm = (first, last, sa, checkpoints)

    #j = getRowBySA(fm, index-1, a, alphabet)
    if index > 0:
        j = sa.index(index-1)

        j2 = LF(fm, newRow)
        while not j == j2:
            newJ = LF(fm, j)
            fm = moveRow(fm, j, j2)
        
            j = newJ
            j2 = LF(fm, j2)

    return fm 

def delete(fm, index):
    ''' Update the BWT for a deletion of a character from position index in the original string '''
    (first, last, sa, checkpoints) = fm

    # Update old row
    row = sa.index(index+1)
    delRow = LF(fm, row)
    remC = last[row]
    tempC = last[delRow]
    last[row] = tempC

    # Delete extra row
    last = last[:delRow] + last[delRow+1:]

    # Update sa
    sa = sa[:delRow] + sa[delRow+1:]
    for i in xrange(len(sa)):
        if sa[i] >= index:
            sa[i] -= 1

    # Update first column
    for k in alphabet:
        if k == remC:
            first[k] = (first[k][0], first[k][1]-1)
        elif k > remC:
            first[k] = (first[k][0]-1, first[k][1]-1)

    # Remove last checkpoint row
    if (len(last)+1) % b == 0:
        checkpoints = checkpoints[:-1]

    # Update checkpoints        
    indexMoved = alphabet.index(tempC)
    indexRemoved = alphabet.index(remC)    
    for x in xrange(row + b - (row%b+1), len(last), b):
        checkpoints[(x+1) / b - 1][indexRemoved] -= 1
    for x in xrange(delRow + b - (delRow%b+1), len(last)-1, b):
        checkpoints[(x+1) / b - 1][alphabet.index(last[x+1])] += 1
    if row < delRow:
        for x in xrange(row + b - (row%b+1), delRow, b):
            checkpoints[(x+1) / b - 1][indexMoved] += 1
    else:
        for x in xrange(delRow + b - (delRow%b+1), row, b):
            checkpoints[(x+1) / b - 1][indexMoved] -= 1


    fm = (first, last, sa, checkpoints)

    j = sa.index(index-1)
    j2 = LF(fm, row)
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
    return (first[c][0] + getCount(last, index, checkpoints, b)) % len(last)

def constructFM(t, b):
    last, sa = bwtViaBwm(t)
    _, tots = rankBwt(last)
    first = firstCol(tots)
    checkpoints = getCheckpoints(last, b)
    return (first, last, sa, checkpoints)
   


for n in xrange(1000):
    print 'Test ' + str(n+1)
    length = random.randint(30,50)
    t = ['$']
    for i in range(length):
        t = [random.choice(['A', 'C', 'T', 'G'])] + t
    print '  ' + ''.join(t)

    alphabet = ['$', 'A', 'C', 'G', 'T']
    b = 3


    fm = constructFM(t, b)


    #insertId = random.randint(0,length)
    #newChar = random.choice(['A', 'C', 'T', 'G'])
    deleteId = random.randint(0,length)

    #t2 = t[:insertId] + [newChar] + t[insertId:]
    t2 = t[:deleteId] + t[deleteId+1:]
    print '  ' + ''.join(t2)

    fm2 = constructFM(t2,b)

    fm_new = delete(fm, deleteId)

#    for i in xrange(len(fm2[1])):
#        print fm2[1][i] + '\t' + str(fm2[2][i]) + '\t' + fm_new[1][i] + '\t' + str(fm_new[2][i])


    for i in xrange(len(fm2[0])):
        if not fm_new[2][i] == fm2[2][i]:
            print 'Error!'
            print fm2[0]
            print fm2[1]
            print fm2[2]
            print fm2[3]

            print ''
            print fm_new[0]
            print fm_new[1]
            print fm_new[2]
            print fm_new[3]
            
            for x in xrange(len(fm2[0])):
                print fm2[0][x] + '\t' + str(fm2[2][x]) + '\t' + fm_new[0][x] + '\t' + str(fm_new[2][x])
            exit()

print 'All correct!'


'''
t = 'TTAGCCAATGGAATGGAAGCCGAT$'
query1 = 'AGCC'
query2 = 'AAT'
query3 = 'ACTGG'
alphabet = sorted(set(t))

b = 4
last, sa = bwtViaBwm(t)
_, tots = rankBwt(last)
first = firstCol(tots)
checkpoints = getCheckpoints(last, b)
fm = (first, last, sa, checkpoints)

print 'T = ' + t
print 'Searching for: ' + query1
matches = find(fm, b, query1)
print matches

print 'Searching for: ' + query2
matches = find(fm, b, query2)
print matches

print 'Searching for: ' + query3
matches = find(fm, b, query3)
print matches



t2 = 'TTAGCCAACTGGAATGGAAGCCGAT$'

last2, sa2 = bwtViaBwm(t2)
_, tots2 = rankBwt(last2)
first2 = firstCol(tots2)
checkpoints2 = getCheckpoints(last2, b)
fm2 = (first2, last2, sa2, checkpoints2)

print 'T = ' + t2
print 'Searching for: ' + query1
matches = find(fm2, b, query1)
print matches

print 'Searching for: ' + query2
matches = find(fm2, b, query2)
print matches

print 'Searching for: ' + query3
matches = find(fm2, b, query3)
print matches


print '\nUpdating initial index'
fm_updated = insert(fm, 8, 'C')
print 'T = ' + t2
print 'Searching for: ' + query1
matches = find(fm_updated, b, query1)
print matches

print 'Searching for: ' + query2
matches = find(fm_updated, b, query2)
print matches

print 'Searching for: ' + query3
matches = find(fm_updated, b, query3)
print matches


(first, last, sa, checkpoints) = fm_updated
for i in xrange(len(last)):
    for c in alphabet:
        if i >= first[c][0] and i < first[c][1]:
            print c + '\t' + last[i] + '\t' + str(sa[i]) + '\t' + str(getCount(last, i, checkpoints, b))

print '\n'
(first, last, sa, checkpoints) = fm2
for i in xrange(len(last)):
    for c in alphabet:
        if i >= first[c][0] and i < first[c][1]:
            print c + '\t' + last[i] + '\t' + str(sa[i]) + '\t' + str(getCount(last, i, checkpoints, b))
'''
