### Team Members: 
Natalie Martinez, Anand Koshy, Nick Garza, Julia Oppenheim

### Overview:
Indexing large collections of reads can be incredibly costly in terms of time and space. Traditional methods often surrender efficiency in one of these two areas for that in the other [2]. In our project, we plan to combine aspects of existing indexing data structures and algorithms to balance the tradeoff between time and space efficiency.

### Methods Used:
* De Bruijn Graph
* Muthukrishnan Algorithm
* Gk-Arrays


### Objective:
*  Data structure that allows for search queries on a large collection of sequencing reads. This involves accommodating for queries such as whether a sequence occurs, in which datasets, and all the positions it occurs at.


### Dependency Installation Instructions:
Some imported libraries we used for this project were Pympler to measure the size of our data structures and an efficient C library called libdivsufsort to build suffix arrays. To download Pympler, a simple conda installation will suffice using `conda install pympler`. To download libdivsufsort, follow the instructions on their installation page: https://github.com/y-256/libdivsufsort. Its installation requires `cmake`.

### Instructions for using Gk Array (and other data structures):

```python
from FMindex import FmIndex
from hash_table import hash_table
from gkarray import GkArray

#inout files
files = [file1, file2, ..., fileN]
#k-mer length
k = 3

# Creating Data Structure
h_table = hash_table(files)
fm_index = FmIndex(files)
gk_array = GkArray(files, k)


# Getting occurrences of k-mer
# 
# returns in the following form:
# {file1: [occurence1, occurence2, ...], file 

k = "TTG"
hash_occurrences = h_table.find_sequence(k)
gk_occurrences = gk_array.get_reads(k)
fm_occurrences = fm_index.occurrences(k)

start_indices = fm_index.start_indices
file_map = fm_index.file_map
# returns files that occurences are in
files = fm_index.report_files(fm_occurrences, start_indices, file_map)


```

### Benchmarking:
We created benchmarking tests for the hash table, FM Index, and Gk Array in order to compare their performances. The test are run on read collection sizes of 50kb, 1 mb, and 100 mb. There are two ways to run the benchmarking scripts. The first is to run the "test_suite.sh" bash script which just runs all of the test scripts for consecutively. 

```python
chmod u+x test_suite.sh
./test_suite.sh
```

The other way is to run the each of the benchmarking files for each of the data structures individually. Each of the benchmarking scripts are titled in the following way: *{name of data structure}_bench_{collection size}.py* 
<br>
For example:
```python
python3 gk_bench_50kb.py
```
<br>

*Disclaimer*: Some of these scripts take a long time to run on the 100 mb file input mainly because the libraries for measuring the size of the data structures take a long time to compute. Please view the benchmarking files for instructions on how to skip the memory measurement step if the script is taking too long to run.



### Sources:
[1] Olsen N, et. al: Challenging a bioinformatic tool's ability to detect microbial contaminants using in silico whole genome sequencing data. PeerJ 2017, 5: e3729.

[2] Maróti Z, Boldogkői Z, Tombácz D et al: Evaluation of whole exome sequencing as an alternative to BeadChip and whole genome sequencing in human population genetic analysis. BMC Genomics (2018): 19(778)

[3] Dolle D, Liu Z, Cotten M, Simpson J, Iqbal Z, Durbin R, McCarthy S, Keane T: Using reference-free compressed data structures to analyse sequencing reads from thousands of human genomes. bioRxiv 2016, 060186.

[4] Philippe N, Salson M, Lecroq T, et al:  Querying large read collections in main memory: a versatile data structure. BMC Bioinformatics 2011: 12(242).

[5] Cox A, Bauer M, Jakobi T, Roson G: Large-scale compression of genomic sequence databases with the burrows–wheeler transform. Bioinformatics 2012, 28(11):1415–1419.

[6]  Välimäki M, Mäkinen V: Space-efficient algorithms for document retrieval. In Annual Symposium on Combinatorial Pattern Matching 2007, 205-215.

[7] Fulton B: Suffix Arrays. Algorithmic Alley 2013. 
http://algorithmicalley.com/archive/2013/06/30/suffix-arrays.aspx

[8] @benfulton 2019. SuffixArrays.
https://github.com/benfulton/Algorithmic-Alley/tree/master/AlgorithmicAlley/SuffixArrays

[9] @kvh 2017. Suffix Trees in Python.
https://github.com/kvh/Python-Suffix-Tree.

