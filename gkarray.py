from suffix_array import SuffixTree
from errorChecks import list_occurrences
from naive_suffix_array import suffix_array_ManberMyers

class GkArray:

    def __init__(self, files, k):
        self.k = k
        self.GkSA = []
        self.GkIFA = None
        self.GkCFPS = None
        self.Cr = "" 
        self.modified_SA = []
        self.file_specs = {} # file number : {file_name, read_length, entries, prev_reads}
        self.starts = [0] #indices of starting positions of files in Cr
        
        
        #For Ukonnen suffix tree to work you have to append $ to end of Cr
        self.concatenate_reads(files)
        #SA = suffix_array_ManberMyers(self.Cr)
        SA = SuffixTree(self.Cr + "$").build_suffix_array()
        #SA = suffix_tree.build_suffix_array()

        self.construct_GkSA(SA)
        self.construct_GkIFA_GkCFA()
        self.construct_GkCFPS()

    def concatenate_reads(self, files):
        file_counter = 0
        total_reads = 0
        for file in files:
            self.Cr, entry = list_occurrences(file, self.Cr)
            self.file_specs[file_counter] = entry
            self.file_specs[file_counter]["prev_read_count"] = total_reads
            self.starts.append(len(self.Cr))

            file_counter += 1
            total_reads += entry["entries"]

    # get the file number corresponding to a certain index for lookup in file_specs
    def get_file(self, x):
        """ Find file from Cr index in O(logk) time """ 
        l = 0
        r = len(self.starts)
        mid = (l + r)//2
        prev = mid + 1
        while l <= r: 
    
            mid = (l + r)//2 
            
            # return if x is in tight range
            if self.starts[mid] <= x and x < self.starts[prev] and mid == prev - 1: 
                return mid 
    
            # If x is greater, ignore left half 
            elif self.starts[mid] < x: 
                l = mid + 1
                prev = mid + 1
    
            # If x is smaller, ignore right half 
            else: 
                r = mid - 1
                prev = mid - 1
        
        # If we reach here, then the element at the edge
        return mid

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

        for i in range (self.k):
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

    # Paper Query 4
    # Given a desired kmer where the position is known, count all occurences of kmer
    def get_total_occurence_count(self, kmer):
        kmer_index = self.find_kmer(kmer)
        if kmer_index < 0:
            return 0
        t = self.GkIFA[kmer_index]
        return self.GkCFPS[t] - self.GkCFPS[t - 1]

    # Paper Query 3
    # Given a desired kmer where the position is known, return a list of indices of occurences
    def get_all_positions(self, kmer):
        kmer_index = self.find_kmer(kmer)
        if kmer_index < 0:
            return []
        t = self.GkIFA[kmer_index]
        upper = self.GkCFPS[t]
        lower = self.GkCFPS[t - 1]
        return self.GkSA[lower : upper + 1]

    #Paper Query 1 and 2
    # Given a desired kmer where the position is known, return a list of reads
    # Answer Query 2 by getting the length of the result returned
    def get_reads(self, kmer):
        l = self.get_all_positions(kmer)
        result = {}
        cardinality = 0
        for i in l:
            file_num = self.get_file(self.g_inverse(i))
            m = self.file_specs[file_num]["read_length"]
            read_num = (self.g_inverse(i) - self.starts[file_num]) // m
            file_name = self.file_specs[file_num]["file"]
            if file_name in result:
                if read_num not in result[file_name]:
                    result[file_name].append(read_num)
                    cardinality += 1
            else:
                result[file_name] = [read_num]
                cardinality += 1
        return result, cardinality

    # Inputs: pattern
    # Returns: index in text if found, -1 if not found (shouldn't happen)
    def find_kmer(self, p):
        t = self.Cr + '$' # t already has terminator
 
        if len(t) == 1: return 1
        l, r = 0, len(self.modified_SA) # invariant: sa[l] < p < sa[r]
        while True:
            c = (l + r) // 2
            # determine whether p < T[sa[c]:] by doing comparisons
            # starting from left-hand sides of p and T[sa[c]:]
            plt = True # assume p < T[sa[c]:] until proven otherwise
            i = 0
            while i < len(p) and self.modified_SA[c]+i < len(t):
                if p[i] < t[self.modified_SA[c]+i]:
                    break # p < T[sa[c]:]
                elif p[i] > t[self.modified_SA[c]+i]:
                    plt = False
                    break # p > T[sa[c]:]
                i += 1 # tied so far
            if plt:
                if c == l + 1:
                    if c < len(t):
                        return self.g(c) # return index transformed to original text
                    else:
                        return -1 # not found, not expected
                r = c
            else:
                if c == r - 1:
                    if c < len(t):
                        return self.g(r) # return index transformed to original text
                    else:
                        return -1 # not found, not expected
                l = c


