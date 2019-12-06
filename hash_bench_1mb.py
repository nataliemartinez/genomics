import glob
from pympler import asizeof
from timeit import default_timer as timer
from hash_table import hash_table

def test_1mb():
    
    h_stats = []
    files = []

    #Collect all file names in collection
    for filename in glob.glob('./input_files/1mb_genome/*.fastq'):
        files.append(filename)

    #building hash table time and memory test

    h_start = timer()
    h_table = hash_table(files)
    h_end = timer()


    h_memory = asizeof.asizeof(h_table)

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

test_1mb()