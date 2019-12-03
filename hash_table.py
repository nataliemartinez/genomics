import sys 
import glob

def __init__(self, files):
    self.reads_dict = {}
    for f in files: 
        self.list_occurrences(f)
    
"""
List occurences of sequences in files
 fh = name of file 
 t = threshold for avg base quality 
 reads_dict = dictionary maping {read: {filename: [line occurrences in file]}}
"""
def list_occurrences(self, f):

    with open(f, 'r') as fh:
        while True:
            first_line = fh.readline()
            if len(first_line) == 0:
                break  # end of file
            name = first_line[1:].rstrip() # line 1 contains name
            seq = fh.readline().rstrip() # line 2 contains sequence 
            if first_pass == True: 
                read_length = len(seq)
                first_pass = False 

            fh.readline()  # line 3 contains '+' so we skip
            qual = fh.readline().rstrip() # line 4 contains quality

            file_occ_dict = self.reads_dict.get(seq, {})
            file_occ = file_occ_dict.get(fh.name, []) 
            file_occ.append(name)
            file_occ_dict[fh.name] = file_occ
            self.reads_dict[seq] = file_occ_dict
    
def find_sequence(self, seq): 
    occurrences = [] #mapping of file to list of reads that k-mer occurs in
    for key in self.reads_dict.keys(): 
        if seq in key: 
            occurrences.append(self.reads_dict.get(key))
    return occurrences


# find number of occurences from all files
def total_reads(self, seq):
    file_occurrences = find_sequence(self, seq)
    total_occurrences = 0
    for f in file_occurrences:
        total_occurrences += len(f.values())
    return total_occurrences