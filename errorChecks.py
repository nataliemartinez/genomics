import sys 

def phred33_to_q(qual):
    return ord(qual)-33

# fh = name of file 
# t = threshold for avg base quality 
# reads_dict = dictionary maping {read: {filename: [line occurrences in file]}}
def listOccurrences(fh, t, reads_dict):

    while True:
        first_line = fh.readline()
        if len(first_line) == 0:
            break  # end of file
        name = first_line[1:].rstrip()
        seq = fh.readline().rstrip()
        fh.readline()  # ignore line starting with +
        qual = fh.readline().rstrip()

        qual_sum = 0
        for i in len(qual): 
            qual_sum += phred33_to_q(qual[i])
        if qual_sum/len(qual) >= t: 
            file_occ_dict = reads_dict.get(seq, {})
            file_occ = file_occ_dict.get(fh, []) 
            file_occ.append(name)
            file_occ_dict[fh] = file_occ
            reads_dict[seq] = file_occ_dict
    
    return reads_dict
    

# def readFiles(): 
#     k_mer_counts = {} 
#     files = sys.argv[1:]
#     for f in files: 
#         if f.find(".fastq") == -1: 
#             print("TODO")
#             # 