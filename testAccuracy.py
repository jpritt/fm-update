#! /usr/bin/env python
import bwt
import random
import time

''' Compare the time to rebuild the index vs updating for substitute, insert, and delete '''

t = 'ACGTAGCGCGTGA$'
alphabet = ['$', 'A', 'C', 'G', 'T']
b = 5

fm = bwt.constructFM(t, b, alphabet)

query = 'CGTAG'

print t
print query

print '0 edits:'
print bwt.findApproximate(fm, b, alphabet, query, 0)
print ''

print '1 edit:'
print bwt.findApproximate(fm, b, alphabet, query, 1)
print ''

print '2 edits:'
print bwt.findApproximate(fm, b, alphabet, query, 2)
print ''
