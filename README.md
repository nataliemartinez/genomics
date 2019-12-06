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


### Sources:
[1] Niko Välimäki and Veli Mäkinen. Space-efficient algorithms for document retrieval. In
Annual Symposium on Combinatorial Pattern Matching, pages 205–215. Springer, 2007.

[2] Anthony J Cox, Markus J Bauer, Tobias Jakobi, and Giovanna Rosone. Large-scale compression of genomic sequence databases with the burrows–wheeler transform. Bioinformatics, 28(11):1415–1419, 2012.

[3] Philippe Nicolas, Salson Mikaël, Lecroq Thierry, Léonard Martine, Commes Thérѐse, and Rivals Eric. Querying large read collections in main memory: a versatile data structure. BMC Bioinformatics, 12.

[4] Dirk-Dominic Dolle, Zhicheng Liu, Matthew L Cotten, Jared T Simpson, Zamin Iqbal, Richard Durbin, Shane McCarthy, and Thomas Keane. Using reference-free compressed data structures to analyse sequencing reads from thousands of human genomes. bioRxiv, page 060186, 2016.

