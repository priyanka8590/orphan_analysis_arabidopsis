import os, sys
import itertools
from itertools import combinations, permutations

def countWaysToSplit(s):
    #split_strings = [''.join(l) for i in range(len(s)) for l in combinations(s, i+1)]
    #split_strings = [''.join(l) for i in range(len(s)) for l in combinations(s, i+1)]
    i = 0
    while i < (len(s)-2):
        a = s[0:i+1]
        i+=1
        print (a)
    j = 0
    while j < (len(s)-1):
        b = s[1:j+1]
        c = s[j+1:len(s)]
        j+=1
        print (b, c)
        
        
    
    """a = s[0]
    j = 1
    while j < (len(s)-1):
        b = s[1:j+1]
        c = s[j+1:len(s)]
        j+=1
        print (a, b, c)
    """
    
countWaysToSplit('xzxzx')