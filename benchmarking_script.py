import glob
import tracemalloc
from timeit import default_timer as timer
from gkarray import GkArray
from FMindex import FmIndex
from hash_table import hash_table
#import unittest


def test_100kb():

    h_stats = []
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


    # #building FM Index time and memory test
    # tracemalloc.start()
    # fm_start = timer()
    # fm_index = FmIndex(files)
    # start_indices = fm_index.start_indices
    # file_map = fm_index.file_map
    # fm_memory = tracemalloc.get_tracemalloc_memory()
    # fm_end = timer()
    # tracemalloc.stop()

    # #query time
    # q_start = timer()
    # occs = fm_index.occurrences("TTG")
    # files = fm_index.report_files(occs, start_indices, file_map)
    # q_end = timer()


    # fm_time = fm_end - fm_start
    # q_time = q_end - q_start
    # fm_stats.extend(((str(fm_time) + " sec"), (str(fm_memory / 1000) + " kb"), (str(q_time) + " sec")))


    # #building Gk Array time and memory test
    # tracemalloc.start()
    # gk_start = timer()
    # gk_array = GkArray(files, 3)
    # gk_memory = tracemalloc.get_tracemalloc_memory()
    # gk_end = timer()
    # tracemalloc.stop()

    # #query time
    # q_start = timer()
    # gk_array.get_reads("TTG")
    # q_end = timer()

    # gk_time = gk_end - gk_start
    # q_time = q_end - q_start
    # gk_stats.extend(((str(gk_time) + " sec"), (str(gk_memory / 1000) + " kb"), (str(q_time) + " sec")))



    # Print performance stats

    titles = ['', 'Hash Table']
    names = ['Build Time', 'Memory', 'Query Time']
    data = [titles] + list(zip(names, h_stats))

    print ("50 kb Collection Results:")
    for i, d in enumerate(data):
        line = '|'.join(str(x).ljust(30) for x in d)
        print(line)
        if i == 0:
            print('-' * len(line))

def test_1mb():
    
    h_stats = []
    files = []

    #Collect all file names in collection
    for filename in glob.glob('./input_files/1mb_genome/*.fastq'):
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

    titles = ['', 'Hash Table']
    names = ['Build Time', 'Memory', 'Query Time']
    data = [titles] + list(zip(names, h_stats))

    print ("1 MB Collection Results:")
    for i, d in enumerate(data):
        line = '|'.join(str(x).ljust(30) for x in d)
        print(line)
        if i == 0:
            print('-' * len(line))


def test_100mb():
    h_stats = []
    files = []

    #Collect all file names in collection
    for filename in glob.glob('./input_files/98mb_genome/*.fastq'):
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

    titles = ['', 'Hash Table']
    names = ['Build Time', 'Memory', 'Query Time']
    data = [titles] + list(zip(names, h_stats))

    print ("100 MB Collection Results:")
    for i, d in enumerate(data):
        line = '|'.join(str(x).ljust(30) for x in d)
        print(line)
        if i == 0:
            print('-' * len(line))


test_100kb()
print ('\n\n')
test_1mb()
print ('\n\n')
test_100mb()
