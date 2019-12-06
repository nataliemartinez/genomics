import sys
import glob
import ctypes
from timeit import default_timer as timer
from suffix_array import SuffixTree
from divsufsort import SuffixArray, divsufsort, divbwt

def CsuffixArray(s):
    """ C library version 10x fast """
    return divsufsort(s, (ctypes.c_int * len(s))())

def CbwtFromSa(t, sa=None):
    """ C library """
    return divbwt(t)

def suffixArray(s):
    """ Our version """
    return SuffixTree(s + '$').build_suffix_array()

def bwtFromSa(t, sa=None):
    ''' Given T, returns BWT(T) by way of the suffix array. '''
    bw = []
    dollarRow = None
    if sa is None:
        sa = suffixArray(t)
    for si in sa:
        if si == 0:
            dollarRow = len(bw)
            bw.append('$')
        else:
            bw.append(t[si-1])
    return ''.join(bw), dollarRow # return string-ized version of list bw

class FmCheckpoints(object):
    ''' Manages rank checkpoints and handles rank queries, which are
        O(1) time, with the checkpoints taking O(m) space, where m is
        length of text. '''
    
    def __init__(self, bw, cpIval=4):
        ''' Scan BWT, creating periodic checkpoints as we go '''
        self.cps = {}        # checkpoints
        self.cpIval = cpIval # spacing between checkpoints
        tally = {}           # tally so far
        # Create an entry in tally dictionary and checkpoint map for
        # each distinct character in text
        for c in bw:
            if c not in tally:
                tally[c] = 0
                self.cps[c] = []
        # Now build the checkpoints
        for i, c in enumerate(bw):
            tally[c] += 1 # up to *and including*
            if i % cpIval == 0:
                for c in tally.keys():
                    self.cps[c].append(tally[c])
    
    def rank(self, bw, c, row):
        ''' Return # c's there are in bw up to and including row '''
        if row < 0 or c not in self.cps:
            return 0
        i, nocc = row, 0
        # Always walk to left (up) when calculating rank
        while (i % self.cpIval) != 0:
            if bw[i] == c:
                nocc += 1
            i -= 1
        return self.cps[c][i // self.cpIval] + nocc

class FmIndex():
    ''' O(m) size FM Index, where checkpoints and suffix array samples are
        spaced O(1) elements apart.  Queries like count() and range() are
        O(n) where n is the length of the query.  Finding all k
        occurrences of a length-n query string takes O(n + k) time.
        
        Note: The spacings in the suffix array sample and checkpoints can
        be chosen differently to achieve different bounds. '''
    
    @staticmethod
    def downsampleSuffixArray(sa, n=4):
        ''' Take only the suffix-array entries for every nth suffix.  Keep
            suffixes at offsets 0, n, 2n, etc with respect to the text.
            Return map from the rows to their suffix-array values. '''
        ssa = {}
        for i, suf in enumerate(sa):
            # We could use i % n instead of sa[i] % n, but we lose the
            # constant-time guarantee for resolutions
            if suf % n == 0:
                ssa[i] = suf
        return ssa
    
    def __init__(self, file_list, cpIval=4, ssaIval=4):
        
        # concat_reads = []

        # for filename in file_list:
        # # do stuff
        #     seq = self.readFile(filename)
        #     self.concat_reads.append(seq)

        self.concat_reads = []
        self.start_indices = []
        self.file_map = {}
        self.start = 0
        self.file_counter = 0
        for filename in file_list:
        # do stuff
            seq = self.readFile(filename)
            self.concat_reads.append(seq)
            
            self.start_indices.append(self.start)
            self.start += len(seq)

            self.file_map[self.file_counter] = filename
            self.file_counter += 1

        self.start_indices.append(self.start) # add end of file index
        self.Cr = ''.join(self.concat_reads)


        if self.Cr[-1] != '$':
            self.Cr += '$' # add dollar if not there already
        # Get BWT string and offset of $ within it
        sa = self.suffixArray(self.Cr)
        self.bwt, self.dollarRow = self.bwtFromSa(self.Cr, sa)
        # Get downsampled suffix array, taking every 1 out of 'ssaIval'
        # elements w/r/t T
        self.ssa = self.downsampleSuffixArray(sa, ssaIval)
        self.slen = len(self.bwt)
        # Make rank checkpoints
        self.cps = FmCheckpoints(self.bwt, cpIval)
        # Calculate # occurrences of each character
        tots = dict()
        for c in self.bwt:
            tots[c] = tots.get(c, 0) + 1
        # Calculate concise representation of first column
        self.first = {}
        totc = 0
        for c, count in sorted(tots.items()):
            self.first[c] = totc
            totc += count

        #return self.start_indices, self.file_map

    def suffixArray(self, s):
        """ Our version """
        return SuffixTree(s + '$').build_suffix_array()

    def bwtFromSa(self, t, sa=None):
        ''' Given T, returns BWT(T) by way of the suffix array. '''
        bw = []
        dollarRow = None
        if sa is None:
            sa = self.suffixArray(t)
        for si in sa:
            if si == 0:
                dollarRow = len(bw)
                bw.append('$')
            else:
                bw.append(t[si-1])
        return ''.join(bw), dollarRow # return string-ized version of list bw
    
    def count(self, c):
        ''' Return number of occurrences of characters < c '''
        if c not in self.first:
            # (Unusual) case where c does not occur in text
            for cc in sorted(self.first.items()):
                if c < cc: return self.first[cc]
            return self.first[cc]
        else:
            return self.first[c]
    
    def range(self, p):
        ''' Return range of BWM rows having p as a prefix '''
        l, r = 0, self.slen - 1 # closed (inclusive) interval
        for i in range(len(p)-1, -1, -1): # from right to left
            l = self.cps.rank(self.bwt, p[i], l-1) + self.count(p[i])
            r = self.cps.rank(self.bwt, p[i], r)   + self.count(p[i]) - 1
            if r < l:
                break
        return l, r+1
    
    def resolve(self, row):
        ''' Given BWM row, return its offset w/r/t T '''
        def stepLeft(row):
            ''' Step left according to character in given BWT row '''
            c = self.bwt[row]
            return self.cps.rank(self.bwt, c, row-1) + self.count(c)
        nsteps = 0
        while row not in self.ssa:
            row = stepLeft(row)
            nsteps += 1
        return self.ssa[row] + nsteps
    
    def hasSubstring(self, p):
        ''' Return true if and only if p is substring of indexed text '''
        l, r = self.range(p)
        return r > l
    
    def hasSuffix(self, p):
        ''' Return true if and only if p is suffix of indexed text '''
        l, r = self.range(p)
        off = self.resolve(l)
        return r > l and off + len(p) == self.slen-1
    
    def occurrences(self, p):
        ''' Return offsets for all occurrences of p, in no particular order '''
        l, r = self.range(p)
        return [ self.resolve(x) for x in range(l, r) ]

    ####### Augmented section ########
    
    def report_files(self, occs, si, file_map):
        """ Report files that occs are from """
        files = []
        for o in occs:
            files.append(file_map[self.find_file_binSearch(si, o)])

        return files

    # def find_file(self, si, x): 
    #     """ find file from Cr index in O(k) time """
    #     for i in range(len(si) - 1):
    #         if x >= si[i] and x < si[i + 1]:
    #             return i
    #     return -1
    
    def find_file_binSearch(self, arr, x):
        """ Find file from Cr index in O(logk) time """ 
        l = 0
        r = len(arr)
        mid = (l + r)//2
        prev = mid + 1
        while l <= r: 
    
            mid = (l + r)//2 
            
            # return if x is in tight range
            if arr[mid] <= x and x < arr[prev] and mid == prev - 1: 
                return mid 
    
            # If x is greater, ignore left half 
            elif arr[mid] < x: 
                l = mid + 1
                prev = mid + 1
    
            # If x is smaller, ignore right half 
            else: 
                r = mid - 1
                prev = mid - 1
        
        # If we reach here, then the element at the edge
        return mid

    
    def readFile(self, f):
        """ Read a fastq in and return its concatenated reads """
        sequence = []
        with open(f, 'r') as fh:
            while True:
                first_line = fh.readline()
                if len(first_line) == 0:
                    break  # end of file
                _ = first_line[1:].rstrip() # line 1 contains name
                seq = fh.readline().rstrip() # line 2 contains sequence
                fh.readline()  # line 3 contains '+' so we skip
                _ = fh.readline().rstrip() # line 4 contains quality

                sequence.append(seq)
        
        return ''.join(sequence)

    
    def buildIndex(self, file_list):
        """ This function grabs all input fastq files, returns FMIndex, start_indices, and file index """ 
        concat_reads = []
        start_indices = []
        file_map = {}
        start = 0
        file_counter = 0
        for filename in file_list:
        # do stuff
            seq = self.readFile(filename)
            concat_reads.append(seq)
            
            start_indices.append(start)
            start += len(seq)

            file_map[file_counter] = filename
            file_counter += 1

        start_indices.append(start) # add end of file index
    

        Cr = ''.join(concat_reads)
        fm = FmIndex(Cr)

        return fm, start_indices, file_map

# def benchmark():
#     """ Rough benchmarking method """
#     start = timer()
#     # Needs to pass in $ to call method from it
#     FM, start_indices, file_map = FmIndex('$').buildIndex(glob.glob('./input_files/1mb_genome/*.fastq'))
#     end = timer()
#     print(str(end - start) + " seconds to build index")

#     start = timer()
#     occs = FM.occurrences("ACCCC")
#     files = FM.report_files(occs, start_indices, file_map)
#     end = timer()
#     print(str(end - start) + " seconds to report occurences")

# benchmark()