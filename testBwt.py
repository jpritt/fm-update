#! /usr/bin/env python
import bwt
import random

for n in xrange(1000):
    print 'Test ' + str(n+1)
    
    # Generate a long random string of random length
    length = random.randint(30,50)
    t = ['$']
    for i in range(length):
        t = [random.choice(['A', 'C', 'T', 'G'])] + t
    print '  ' + ''.join(t)

    alphabet = ['$', 'A', 'C', 'G', 'T']
    b = 3

    # Construct the fm index
    fm = bwt.constructFM(t, b, alphabet)


    letters = set(t)
    letters.remove('$')
    
    # Test substitution of a character
    #subId = random.randint(0,length-1)
    #newChar = random.choice(list(letters))
    #t2 = t[:subId] + [newChar] + t[subId+1:]

    # Test insertion of a character
    #insertId = random.randint(0,length-1)
    #newChar = random.choice(['A', 'C', 'T', 'G'])
    #t2 = t[:insertId] + [newChar] + t[insertId:]

    # Test deletion of a character
    deleteId = random.randint(0,length-1)
    t2 = t[:deleteId] + t[deleteId+1:]
    
    print '  ' + ''.join(t2)

    fm2 = bwt.constructFM(t2, b, alphabet)

    #fm_new = bwt.substitute(fm, b, alphabet, subId, newChar)
    #fm_new = bwt.insert(fm, b, alphabet, insertId, newChar)
    fm_new = bwt.delete(fm, b, alphabet, deleteId)


    # Check results for accuracy
    for i in xrange(len(fm2[0])):
        if not fm_new[2][i] == fm2[2][i] and fm_new[1][i] == fm2[1][i]:
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

