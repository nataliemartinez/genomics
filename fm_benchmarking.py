import glob
from timeit import default_timer as timer
from pympler import asizeof
from FMindex import FmIndex
#import unittest


def test_100kb():

    fm_stats = []
    files = []

    #Collect all file names in collection
    for filename in glob.glob('./input_files/50kb_genome/*.fastq'):
        files.append(filename)


    # #building FM Index time and memory test
    fm_start = timer()
    fm_index = FmIndex(files)
    start_indices = fm_index.start_indices
    file_map = fm_index.file_map
    fm_memory = asizeof.asizeof(fm_index)
    fm_end = timer()


    #query time
    q_start = timer()
    occs = fm_index.occurrences("TTG")
    files = fm_index.report_files(occs, start_indices, file_map)
    q_end = timer()


    fm_time = fm_end - fm_start
    q_time = q_end - q_start
    fm_stats.extend(((str(fm_time) + " sec"), (str(fm_memory / 1000) + " kb"), (str(q_time) + " sec")))



    # Print performance stats

    titles = ['', 'FM Index']
    names = ['Build Time', 'Memory', 'Query Time']
    data = [titles] + list(zip(names, fm_stats))

    print ("50 kb Collection Results:")
    for i, d in enumerate(data):
        line = '|'.join(str(x).ljust(30) for x in d)
        print(line)
        if i == 0:
            print('-' * len(line))

def test_1mb():
    
    fm_stats = []
    files = []

    #Collect all file names in collection
    for filename in glob.glob('./input_files/1mb_genome/*.fastq'):
        files.append(filename)


    # #building FM Index time and memory test
    fm_start = timer()
    fm_index = FmIndex(files)
    start_indices = fm_index.start_indices
    file_map = fm_index.file_map
    fm_memory = asizeof.asizeof(fm_index)
    fm_end = timer()


    #query time
    q_start = timer()
    occs = fm_index.occurrences("TTG")
    files = fm_index.report_files(occs, start_indices, file_map)
    q_end = timer()


    fm_time = fm_end - fm_start
    q_time = q_end - q_start
    fm_stats.extend(((str(fm_time) + " sec"), (str(fm_memory / 1000) + " kb"), (str(q_time) + " sec")))


    # Print performance stats

    titles = ['', 'FM Index']
    names = ['Build Time', 'Memory', 'Query Time']
    data = [titles] + list(zip(names, fm_stats))

    print ("1 MB Collection Results:")
    for i, d in enumerate(data):
        line = '|'.join(str(x).ljust(30) for x in d)
        print(line)
        if i == 0:
            print('-' * len(line))


def test_100mb():
    fm_stats = []
    files = []

    #Collect all file names in collection
    for filename in glob.glob('./input_files/98mb_genome/*.fastq'):
        files.append(filename)


    #building Gk Array time and memory test
    # #building FM Index time and memory test
    fm_start = timer()
    fm_index = FmIndex(files)
    #start_indices = fm_index.start_indices
    #file_map = fm_index.file_map
    fm_memory = asizeof.asizeof(fm_index)
    fm_end = timer()


    #query time
    q_start = timer()
    #occs = fm_index.occurrences("TTG")
    #files = fm_index.report_files(occs, start_indices, file_map)
    q_end = timer()


    fm_time = fm_end - fm_start
    q_time = q_end - q_start
    fm_stats.extend(((str(fm_time) + " sec"), (str(fm_memory / 1000) + " kb"), (str(q_time) + " sec")))


    # Print performance stats

    titles = ['', 'FM Index']
    names = ['Build Time', 'Memory', 'Query Time']
    data = [titles] + list(zip(names, fm_stats))

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