# Suffix Array Implementation that runs in O(nlogn) 
# using Manber-Myers Algorithm

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

def suffix_array_ManberMyers(str):
    return sort_bucket(str, (i for i in range(len(str))))

if __name__ == "__main__":
    with open("input_files/input.txt") as f:
        m = f.read()
    str = "aacaactcaattcaaacaag"
    output = suffix_array_ManberMyers(str)
    print (output)