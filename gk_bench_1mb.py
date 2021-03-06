import glob
from timeit import default_timer as timer
from gkarray import GkArray
from pympler import asizeof

def test_1mb():
    
    gk_stats = []
    files = []

    #Collect all file names in collection
    for filename in glob.glob('./input_files/1mb_genome/*.fastq'):
        files.append(filename)


    #building Gk Array time and memory test

    gk_start = timer()
    gk_array = GkArray(files, 3)
    gk_end = timer()


    gk_memory = asizeof.asizeof(gk_array)
    
    #query time
    q_start = timer()
    gk_array.get_reads("TTG")
    q_end = timer()

    gk_time = gk_end - gk_start
    q_time = q_end - q_start
    gk_stats.extend(((str(gk_time) + " sec"), (str(gk_memory / 1000) + " kb"), (str(q_time) + " sec")))



    # Print performance stats

    titles = ['', 'Gk Array']
    names = ['Build Time', 'Memory', 'Query Time']
    data = [titles] + list(zip(names, gk_stats))

    print ("1 MB Collection Results:")
    for i, d in enumerate(data):
        line = '|'.join(str(x).ljust(30) for x in d)
        print(line)
        if i == 0:
            print('-' * len(line))

test_1mb()