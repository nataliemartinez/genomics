class GkArray:

    def __init__(self, reads, k, error_map, read_length, Cr, SA):
        self.reads = reads
        self.k = k
        self.error_map = error_map
        self.read_length = read_length
        self.Cr = None
        self.GkSA = None
        self.GkIFA = None
        self.GkCFPS = None
        self.Cr = Cr

        self.construct_GkSA(SA)
        self.construct_GkIFA_GkCFA()
        self.construct_GkCFPS()
        

    def concatenate_reads(self):
        with open(self.reads) as read_file:
            input_lines = read_file.readlines()
        
        total_reads = int(len(input_lines) / 4)
        Cr = ""

        for i in range(total_reads):
            read = input_lines[i*4][1:].rstrip()
            seq = input_lines[(i*4)+1].rstrip().upper()
            if read in self.error_map:
                Cr += seq

        self.Cr = Cr
    
    def g(self, index):
        m_hat = self.read_length - self.k + 1
        return index // self.read_length * m_hat + (index % self.read_length)

    def is_P_position(self, index):
        return (index % self.read_length) <= (self.read_length - self.k)
    
    def construct_GkSA(self, SA):

        GkSA = []

        for j in range(len(SA)):
            if self.is_P_position(SA[j]):
                GkSA.append(self.g(SA[j]))
                print (str(SA[j]) + ":" + str(self.g(SA[j])))
        
        self.GkSA = GkSA

    def f_q(self, index):
        return index #fix

    def construct_GkIFA_GkCFA(self):
        GkIFA = GkCFA = [0] * len(self.GkSA)
        GkCFA[0] = 1
        t = 0

        for i in range(1, len(self.GkSA)):
            j = self.GkSA[i]
            j_hat = self.GkSA[i - 1]
            if self.f_q(j) != self.f_q(j_hat):
                t += 1
                GkCFA[t] = 0
            GkIFA[j] = t
            GkCFA[t] += 1

        self.GkIFA = GkIFA
        self.GkCFPS = GkCFA #will update in place
    
    def construct_GkCFPS(self):
        for i in range(1, len(self.GkCFPS)):
            self.GkCFPS[i] += self.GkCFPS[i - 1]


def main():

    reads = ["aacaact", "caattca", "aacaagc"]
    SA = [13, 0, 14, 3, 17, 8, 1, 15,4, 18, 9, 20, 12, 2, 16, 7, 5, 19, 11, 6, 10]
    em = {}

    tester = GkArray(reads, 3, em, 7, "aacaactcaattcaaacaagc", SA)

    print (tester.GkSA)
    print (tester.GkIFA)
    print (tester.GkCFPS)


main()




