import glob
import tracemalloc
from timeit import default_timer as timer
from gkarray import GkArray
from hash_table import hash_table
#import unittest


def test_100kb():

    h_stats, gk_stats, fm_stats = ([] for i in range(3))
    files = []

    #Collect all file names in collection
    for filename in glob.glob('./input_files/50kb_genome/*.fastq'):
        files.append(filename)

    #building hash table time and memory test
    tracemalloc.start()
    h_start = timer()
    h_table = hash_table(files)
    h_memory = tracemalloc.get_tracemalloc_memory()
    h_end = timer()
    tracemalloc.stop()

    #query time
    q_start = timer()
    h_table.find_sequence("TTG")
    q_end = timer()


    h_time = h_end - h_start
    q_time = q_end - q_start
    h_stats.extend(((str(h_time) + " sec"), (str(h_memory / 1000) + " kb"), (str(q_time) + " sec")))


    #building FM Index time and memory test
    tracemalloc.start()
    fm_start = timer()
    #fm_index = FmIndex("Build")
    fm_memory = tracemalloc.get_tracemalloc_memory()
    fm_end = timer()
    tracemalloc.stop()

    #query time
    q_start = timer()
    #FM index query
    q_end = timer()


    fm_time = fm_end - fm_start
    q_time = q_end - q_start
    fm_stats.extend(((str(fm_time) + " sec"), (str(fm_memory / 1000) + " kb"), (str(q_time) + " sec")))


    #building Gk Array time and memory test
    tracemalloc.start()
    gk_start = timer()
    gk_array = GkArray(files, 3)
    gk_memory = tracemalloc.get_tracemalloc_memory()
    gk_end = timer()
    tracemalloc.stop()

    #query time
    q_start = timer()
    gk_array.get_reads("TTG")
    q_end = timer()

    gk_time = gk_end - gk_start
    q_time = q_end - q_start
    gk_stats.extend(((str(gk_time) + " sec"), (str(gk_memory / 1000) + " kb"), (str(q_time) + " sec")))



    # Print performance stats

    titles = ['', 'Hash Table', 'FM Index', 'Gk Array']
    names = ['Build Time', 'Memory', 'Query Time']
    data = [titles] + list(zip(names, h_stats, fm_stats, gk_stats))

    print ("50 kb Collection Results:")
    for i, d in enumerate(data):
        line = '|'.join(str(x).ljust(30) for x in d)
        print(line)
        if i == 0:
            print('-' * len(line))

def test_1mb():
    print ("fill in")
    return 0

def test_100mb():
    print ("fill in")
    return 0

# def test_hash():
#     files = []

#     #Collect all file names in collection
#     for filename in glob.glob('./input_files/50kb_genome/*.fastq'):
#         files.append(filename)

#     h_table = hash_table(['input_files/test.txt'])
#     out = h_table.find_sequence("TTG")
#     out1 = h_table.total_reads("TTG")
#     print (out, out1)

# if __name__ == '__main__':
    # print ("out")
    # test_hash_table()
#test_hash()
test_100kb()
print ('\n\n\n')
test_1mb()
print ('\n\n\n')
test_100mb()
