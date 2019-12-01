from suffix_array import SuffixTree
from errorChecks import list_occurrences

class GkArray:

    def __init__(self, files, k):
        self.k = k
        self.valid_reads = {}
        self.GkSA = []
        self.GkIFA = None
        self.GkCFPS = None
        self.Cr = "" 
        self.modified_SA = []
        self.file_specs = {} # file number : {file_name, read_length, entries, prev_reads}
        self.starts = [0] #indices of starting positions of files in Cr

        #For Ukonnen suffix tree to work you have to append $ to end of Cr
        suffix_tree = SuffixTree(self.Cr + "$")
        SA = suffix_tree.build_suffix_array()

        self.concatenate_reads(files)
        self.construct_GkSA(SA)
        self.construct_GkIFA_GkCFA()
        self.construct_GkCFPS()

    def concatenate_reads(self, files):
        file_counter = 0
        total_reads = 0
        for file in files:
            self.valid_reads, self.Cr, entry = list_occurrences(file, self.valid_reads, self.Cr)
            self.file_specs[file_counter] = entry
            self.file_specs[file_counter]["prev_read_count"] = total_reads
            self.starts.append(len(self.Cr))

            file_counter += 1
            total_reads += entry["entries"]
    
    #O(n), could be faster with modified binary search?
    # get the file number corresponding to a certain index for lookup in file_specs
    def get_file(self, index):
        for i in range(len(self.starts) - 1):
            if index >= self.starts[i] and index < self.starts[i + 1]:
                return i
        return None

    # transform index of Cr to index for GkSA, modified for multiple files
    def g(self, index):
        file_num = self.get_file(index)
        entry = self.file_specs[file_num]
        read_num = (index - self.starts[file_num]) // entry["read_length"]
        total_reads_before = entry["prev_read_count"] + read_num
        return index - (self.k - 1) * total_reads_before

    def g_inverse(self, index_in_gksa):
        return self.modified_SA[index_in_gksa]

    # compares two kmers to each other
    def compare_kmer(self, index1, index2):
        g_index1 = self.g_inverse(index1)
        g_index2 = self.g_inverse(index2)

        for i in range (self.k + 1):
            if self.Cr[g_index1 + i] != self.Cr[g_index2 + i]:
                return False
        return True

    # determines whether an index is a valid kmer offset within a read
    def is_P_position(self, index):
        file_num = self.get_file(index)
        read_length = self.file_specs[file_num]["read_length"]
        return ((index - self.starts[file_num]) % read_length) <= (read_length - self.k)
    
    # use suffix array to construct GkSA, per paper specifications
    def construct_GkSA(self, SA):
        for j in range(len(SA)):
            if self.is_P_position(SA[j]):
                self.GkSA.append(self.g(SA[j]))
                self.modified_SA.append(SA[j])

    # use GkSA to construct GkIFA and GkCFA, per paper specifications
    def construct_GkIFA_GkCFA(self):
        GkIFA = [0] * len(self.GkSA)
        GkCFA = [0] * len(self.GkSA)
        GkCFA[0] = 1
        t = 0

        for i in range(1, len(self.GkSA)):
            if not self.compare_kmer(i, i - 1):
                t += 1
                GkCFA[t] = 0
            GkIFA[self.GkSA[i]] = t
            GkCFA[t] += 1
        
        self.GkIFA = GkIFA
        self.GkCFPS = GkCFA[:t + 1] #will update in place

    # use GkCFA to construct GkCFPS, per paper specifications
    def construct_GkCFPS(self):
        for i in range(1, len(self.GkCFPS)):
            self.GkCFPS[i] += self.GkCFPS[i - 1]


    def get_rank(self, read_num, read_index, file_num):
        j = self.g(read_num * self.read_length + read_index + self.starts[file_num])
        return self.GkIFA[j]

    # Paper Query 4
    # Given a desired kmer where the position is known, count all occurences of kmer
    def get_total_occurence_count(self, read_num, read_index, file_num):
        t = self.get_rank(read_num, read_index, file_num)
        return self.GkCFPS[t] - self.GkCFPS[t - 1]

    # Paper Query 3
    # Given a desired kmer where the position is known, return a list of indices of occurences
    def get_all_positions(self, read_num, read_index, file_num):
        t = self.get_rank(read_num, read_index, file_num)
        upper = self.GkCFPS[t]
        lower = self.GkCFPS[t - 1]
        return self.GkSA[lower : upper + 1]

    #Paper Query 1 and 2
    # Given a desired kmer where the position is known, return a list of reads
    # Answer Query 2 by getting the length of the result returned
    def get_reads(self, read_num, read_index, file_num):
        l = self.get_all_positions(read_num, read_index, file_num)
        m = 0 #get length of read for that read?
        result = []
        for i in l:
            result.append(self.g_inverse(i) // m)
        return result




