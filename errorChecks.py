import sys 
import glob

""" Convert ascii char to quality score """
def phred33_to_q(qual):
    return ord(qual)-33

def k_mer_quality(kmer): 
    qual_sum = 0
    for i in range(len(kmer)): # get total of quality scores
        qual_sum += phred33_to_q(kmer[i]) 
    return qual_sum / 3 

"""
List occurences of sequences in files
 fh = name of file 
 t = threshold for avg base quality 
"""
def list_occurrences(f, cr):
    min_qual = sys.maxsize
    first_pass = True #used for getting read length for file 
    read_length = 0
    entries = 0

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

            for i in range(len(seq) - 3): #check this 
                # kmer = seq[i:i+3] --> might be helpful for building Gk array 
                kmer_quals = qual[i:i+3]
                avg_qual = k_mer_quality(kmer_quals)
                if avg_qual < min_qual:
                    min_qual = avg_qual

            # if avg quality > threshold add to dict {seq: {file: [occ,...]}}
            if min_qual >= QUALITY_THRESHOLD:  # quality threshold = 1; should this be inclusive or exclusive? 
                cr += seq
                entries += 1
        
        entry = {"file" : f, "read_length" : read_length, "entries" : entries}
    
    return cr, entry
    

QUALITY_THRESHOLD = 1