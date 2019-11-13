import sys 
import glob

""" Convert ascii char to quality score """
def phred33_to_q(qual):
    return ord(qual)-33

"""
List occurences of sequences in files
 fh = name of file 
 t = threshold for avg base quality 
 reads_dict = dictionary maping {read: {filename: [line occurrences in file]}}
"""
def listOccurrences(f, t, reads_dict):

    with open(f, 'r') as fh:
        while True:
            first_line = fh.readline()
            if len(first_line) == 0:
                break  # end of file
            name = first_line[1:].rstrip() # line 1 contains name
            seq = fh.readline().rstrip() # line 2 contains sequence
            fh.readline()  # line 3 contains '+' so we skip
            qual = fh.readline().rstrip() # line 4 contains quality

            qual_sum = 0
            for i in range(len(qual)): # get total of quality scores
                qual_sum += phred33_to_q(qual[i])

            if qual_sum/len(qual) >= t:  # if avg quality > threshold add to dict {seq: {file: [occ,...]}}
                file_occ_dict = reads_dict.get(seq, {})
                file_occ = file_occ_dict.get(fh.name, []) 
                file_occ.append(name)
                file_occ_dict[fh.name] = file_occ
                reads_dict[seq] = file_occ_dict
    
    return reads_dict
    
"""
Read all the fastqs in our directory and build map {seq: {file: [occ,...]}}
"""
def readFiles(): 
    k_mer_counts = {} # map of {seq: {file: [occ,...]}}

    for filename in glob.glob('*.fastq'):
        # do stuff
        listOccurrences(filename, 5, k_mer_counts)

    for k,v in k_mer_counts.items():
        print(k,end=' ')
        print(v)

readFiles()