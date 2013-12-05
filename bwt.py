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
        for x in bw[i:(nextC+1)]:
            if x == bw[i]:
                count -= 1
        return count

def findRange(fm, b, alphabet, substring, start, end):
    ''' Find the location of substring in the fm index between b-values start and end '''
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
                    currId  = getCount(last, i, checkpoints, b, alphabet)
                    if currId < minId:
                        minId = currId
                    if currId > maxId:
                        maxId = currId
            if minId > maxId:
                return []
            else:
                return findRange(fm, b, alphabet, substring[:-1], int(minId), int(maxId)+1)

def find(fm, b, alphabet, substring):
    ''' Find all indexes of the substring in the text represented in bw '''
    first = fm[0]
    (minI, maxI) = first[substring[-1]]
    return findRange(fm, b, alphabet, substring, 0, maxI-minI)

def moveRow(fm, b, alphabet, i, j):
    ''' Moves the row at index i in the bwt to index j, shifting all rows in-between and updates SA and checkpoints '''
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

def insert(fm, b, alphabet, index, c):
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
    bVal = getCount(last, row, checkpoints, b, alphabet)

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

        j2 = LF(fm, b, alphabet, newRow)
        while not j == j2:
            newJ = LF(fm, b, alphabet, j)
            fm = moveRow(fm, b, alphabet, j, j2)
        
            j = newJ
            j2 = LF(fm, b, alphabet, j2)

    return fm 

def delete(fm, b, alphabet, index):
    ''' Update the BWT for a deletion of a character from position index in the original string '''
    (first, last, sa, checkpoints) = fm

    # Update old row
    row = sa.index(index+1)
    delRow = LF(fm, b, alphabet, row)
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
    for k in first.keys():
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
    for x in xrange(delRow + b - (delRow%b+1), len(last), b):
        checkpoints[(x+1) / b - 1][alphabet.index(last[x])] += 1
    if row < delRow:
        for x in xrange(row + b - (row%b+1), delRow, b):
            checkpoints[(x+1) / b - 1][indexMoved] += 1
    else:
        for x in xrange(delRow + b - (delRow%b+1), row, b):
            checkpoints[(x+1) / b - 1][indexMoved] -= 1

    if row > delRow:
        row -= 1


    fm = (first, last, sa, checkpoints)

    if index > 0:
        j = sa.index(index-1)
        j2 = LF(fm, b, alphabet, row)
        while not j == j2:
            newJ = LF(fm, b, alphabet, j)
            fm = moveRow(fm, b, alphabet, j, j2)
        
            j = newJ
            j2 = LF(fm, b, alphabet, j2)

    return fm 

def substitute(fm, b, alphabet, index, c):
    ''' Update the BWT for a substitution of a new character at position index in the original string '''
    (first, last, sa, checkpoints) = fm

    # Update old row
    row = sa.index(index+1)
    remC = last[row]
    nextRow = LF(fm, b, alphabet, row)
    last[row] = c
    
    # Update first column
    for k in first.keys(): 
        if k == remC:
            first[k] = (first[k][0], first[k][1]-1)
        if k == c:
            first[k] = (first[k][0], first[k][1]+1)
        if k > remC and k <= c:
            first[k] = (first[k][0]-1, first[k][1]-1)
        if k > c and k <= remC:
            first[k] = (first[k][0]+1, first[k][1]+1)

    # Update checkpoints        
    indexAdded = alphabet.index(c)
    indexRemoved = alphabet.index(remC)    
    for x in xrange(row + b - (row%b+1), len(last), b):
        checkpoints[(x+1) / b - 1][indexAdded] += 1
        checkpoints[(x+1) / b - 1][indexRemoved] -= 1

    # New row index for next row
    newRowPos = LF(fm, b, alphabet, row)
    fm = moveRow((first, last, sa, checkpoints), b, alphabet, nextRow, newRowPos)

    if index > 0:
        j = fm[2].index(index-1)
        j2 = LF(fm, b, alphabet, newRowPos)
        while not j == j2:
            newJ = LF(fm, b, alphabet, j)
            fm = moveRow(fm, b, alphabet, j, j2)
       
            j = newJ
            j2 = LF(fm, b, alphabet, j2)

    return fm 

def LF(fm, b, alphabet, index):
    ''' Step forward one step in the bwt and return the next row '''
    (first, last, sa, checkpoints) = fm

    c = last[index]
    return (first[c][0] + getCount(last, index, checkpoints, b, alphabet)) % len(last)

def constructFM(t, b, alphabet):
    ''' Construct an FM-index for string t, with checkpoints every b rows and the predefined alphabet'''
    last, sa = bwtViaBwm(t)
    _, tots = rankBwt(last)
    first = firstCol(tots)
    checkpoints = getCheckpoints(last, b, alphabet)
    return (first, last, sa, checkpoints)

def findApproximate(fm, b, alphabet, substring, k):
    ''' Find all approximate matches to the fm index by making all possible mutated strings and searching for exact matches '''
    variations = dict()
    variations[substring] = []
    for i in xrange(k):
        # first check if any of the current strings match
        matches = dict()
        for k,v in variations.items():
            for m in find(fm, b, alphabet, k):
                if m in matches:
                    matches[m] += k
                else:
                    matches[m] = v
        if len(matches) > 0:
            return matches          

        newV = dict()
        for v in variations:
            for var,edit in makeEdits(v).items():
                e = variations[v] + [edit]
                if not var in newV:
                    newV[var] = e
                else:
                    newV[var] = min(newV[var], e)
        variations = newV

    matches = dict()
    for v,k in variations.items():
        for m in find(fm, b, alphabet, v):
            if m in matches:
                matches[m] += k
            else:
                matches[m] = k
    return matches          


def makeEdits(t):
    ''' Returns all strings at edit distance 1 from t '''
    variations = dict()
    #variations += [t]
    chars = ['A', 'C', 'G', 'T']    

    # 0 = substitution
    # 1 = insertion
    # 2 = deletion

    # deletions
    for i in xrange(len(t)):
        s = t[:i] + t[i+1:]
        if not s in variations:
            variations[s] = (1, i, t[i])
        else:
            variations[s] = min(variations[s], (1, i, t[i]))

    # insertions
    for i in xrange(len(t)+1):
        for c in chars:
            s = t[:i] + c + t[i:]
            if not s in variations:
                variations[s] = (2, i)
            else:
                variations[s] = min(variations[s], (2, i))

    # substitutions
    for i in xrange(len(t)):
        for c in chars:
            if not c == t[i]:
                s = t[:i] + c + t[i+1:]
                if not s in variations:
                    variations[s] = (0, i, t[i])
                else:
                    variations[s] += min(variations[s], (0, i, t[i]))

    return variations

def editDistance(t1, t2):
    ''' Return the minimal sequence of edits to get from t1 to t2 '''
    m = np.zeros([len(t2)+1, len(t1)+1])
    for i in xrange(len(t1)):
        m[0][i+1] = i+1
    for i in xrange(len(t2)):
        m[i+1][0] = i+1

    for x in xrange(len(t2)):
        for y in xrange(len(t1)):
            match = 0 if t2[x] == t1[y] else 1
            m[x+1][y+1] = min(m[x][y+1]+1, m[x+1][y]+1, m[x][y]+match)

#    print m

    return editSeq(m, t1, t2)

def editSeq(m, t1, t2):
    ''' Return the sequence of edits that corresponds to the matrix m '''
    x = len(m) - 1
    y = len(m[0]) - 1
    edits = []

    while x > 0 and y > 0:
#        print '(' + str(x) + ', ' + str(y) + ')'
        match = 0 if t2[x-1] == t1[y-1] else 1
        if m[x][y] == m[x-1][y-1] + match:
            if not t2[x-1] == t1[y-1]:
                edits = [(0, y-1, t2[x-1])] + edits
            x -= 1
            y -= 1
        elif m[x][y] == m[x-1][y] + 1:
            edits = [(1, x-1, t2[x-1])] + edits
            x -= 1
        else:
            edits = [(2, y-1)] + edits
            y -= 1

    if y > 0:
        while y > 0:
            edits = [(2, y-1)] + edits
            y -= 1
    else:
        while x > 0:
            edits = [(1, x, t2[x-1])] + edits
            x -= 1
    return edits


