import glob
from gkarray import GkArray
#from FMindex import FmIndex
from hash_table import hash_table
#import unittest



def test_hash_table():
    #for filename in glob.glob('./input_files/*.fastq'):
    h_table = hash_table(["input_files/test.txt"])
    print (h_table.find_sequence("TGG"))





# if __name__ == '__main__':
    # print ("out")
    # test_hash_table()

test_hash_table()
