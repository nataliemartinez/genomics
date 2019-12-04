import glob
import tracemalloc
from timeit import default_timer as timer
from gkarray import GkArray
#from FMindex import FmIndex
from hash_table import hash_table
#import unittest




def test_100kb():
    #for filename in glob.glob('./input_files/*.fastq'):

    h_stats, gk_stats, fm_stats = ([] for i in range(3))
    files = []
    for filename in glob.glob('./input_files/50kb_genome/*.fastq'):
        files.append(filename)

    #building hash table time and memory test
    tracemalloc.start()
    h_start = timer()
    h_table = hash_table(["input_files/test.txt"])
    h_memory = tracemalloc.get_tracemalloc_memory()
    h_end = timer()
    tracemalloc.stop()

    h_time = h_end - h_start
    h_stats.extend(((str(h_time) + " sec"), (str(h_memory / 1000) + " kb"), 0))


    #building FM Index time and memory test
    tracemalloc.start()
    fm_start = timer()
    fm_index = 0
    fm_memory = tracemalloc.get_tracemalloc_memory()
    fm_end = timer()
    tracemalloc.stop()

    fm_time = fm_end - fm_start
    fm_stats.extend(((str(fm_time) + " sec"), (str(fm_memory / 1000) + " kb"), 0))


    #building Gk Array time and memory test
    tracemalloc.start()
    gk_start = timer()
    gk_array = 0
    gk_memory = tracemalloc.get_tracemalloc_memory()
    gk_end = timer()
    tracemalloc.stop()

    gk_time = gk_end - gk_start
    gk_stats.extend(((str(gk_time) + " sec"), (str(gk_memory / 1000) + " kb"), 0))



    titles = ['','Hash Table', 'FM Index', 'Gk Array']
    names = ['Build Time', 'Memory', 'Query Time']
    data = [titles] + list(zip(names, h_stats, fm_stats, gk_stats))

    for i, d in enumerate(data):
        line = '|'.join(str(x).ljust(30) for x in d)
        print(line)
        if i == 0:
            print('-' * len(line))

def test_1mb():
    return 0

def test_100mb():
    return 0



# if __name__ == '__main__':
    # print ("out")
    # test_hash_table()

test_100kb()
