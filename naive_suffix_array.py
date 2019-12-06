# Suffix Array Implementation that runs in O(nlogn) 
# using Manber-Myers Algorithm
# Source code taken from here:
# https://github.com/benfulton/Algorithmic-Alley/tree/master/AlgorithmicAlley/SuffixArrays

import time
from collections import defaultdict, Counter

def sort_bucket(str, bucket, order=1):
    d = defaultdict(list) 
    for i in bucket:
        key = str[i:i+order]
        d[key].append(i)
    result = []
    for k,v in sorted(d.items()):
        if len(v) > 1:
            result += sort_bucket(str, v, order*2)
        else:
            result.append(v[0])
    return result

def naive_suffix_array(str):
    return sort_bucket(str, (i for i in range(len(str))))
