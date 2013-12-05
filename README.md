fm-update
=========

Final project for Computational Genomics

Uses an FM index to map reads to a reference genome. Continuously updates the FM index to match the genome from which reads were sequenced.

I have included all python files that I used to create and test the algorithms in this project. I changed many of the test scripts as I went to collect data, so a few of them may no longer work. I have listed the most important files, with descriptions of their content, below.

bwt.py
Contains function to create and update the fm-index and search for substrings.
  constructFM - constructs an FM-index for a given string
  insert - update the FM-index for an insertion
  delete - update the FM-index for a deletion
  substitute - update the FM-index for a substitution
  findApproximate - return the indices of all matches for a substring against the fm-index, with up to k errors

iterativeEM.py
Contains one function, iterativeEM, which takes as input an fm-index, a set of reads, and a list of parameters. The function runs the iterative referencing algorithm and returns the accuracy and space of the result. The optional argument prop indicates the proportion of reads to contribute to the list of genome mutations. The default is 1 (all reads contribute), which corresponds to the original iterative referencing algorithm.

iterativeEMDist.py
Contains one function, iterativeEMDist. This function is identical to iterativeEM above, except it stores a coverage vector for mapped reads and allows only the first 'depth' reads that overlap a base to contribute to the mutations list for that base. The parameters depth and chunkSize for this function correspond to the number of contributing reads and the width of the chunks in the coverage vector. Like iterativeEM, this function returns the accuracy and space requirements of the results.

===========================================
The following files were used to test code and generate results 
(approximately in chronological order)
===========================================
testBwt.py
I used this script to make sure that the insert, delete, and substitute functions were correctly updating the fm-index. To test this, I repeatedly generated a long random string and constructed the fm-index for it. I then introduced a mutation and both constructed a new fm-index and updated the old one, and checked to make sure the 2 indexes were identical.

analyzeTime.py
Compares the time to rebuild the fm-index vs. time to update it for varying genome lengths

testAccuracy.py
Computes the accuracy of matching reads to a mutated reference genome for varying parameters. In this case, sequencing errors are not introduced into reads. 

testIterative.py
Compares the accuracy of the iterative referencing algorithm to the initial matching accuracy for varying values of read length, number of reads, and mutation frequency. In this case, sequencing errors are not introduced into reads.

testIterativeError.py
Same as testIterative.py above, but now sequencing errors are introduced into reads. Each nucleotide has some small probability of being replaced by a different base (substitution errors).

testIterativeErrorMargins.py
Same as testIterativeError.py above, but a small margin (default 2*readLen) on each end of the genome is not mutated. This was meant to test whether the algorithm would perform better on genomes with cluster mutations and regions with low-density mutations. However, preliminary testing showed no significant difference.

compareMethods.py
This is a testing script used to generate the comparative results in the paper. The script varies the read length, number of reads, and error frequency and compares the accuracy and space results for the original iterative referencing algorithm, the reduced algorithm (1% of reads stored), and the algorithm with the coverage vector for chunkSize=1--depth=1, chunkSize=1--depth=2, chunkSize=5--depth=2.

compareMethodsWeighted.py
Testing script identical to compareMethods except the generated reads are twice as likely to originate from the first half of the genome. This yields variable coverage to more accurately reflect actual data. The script compares accuracy and space results for the original iterative referencing algorithm, the reduced algorithm (1% of reads stored), and the algorithm with coverage vector for chunkSize=50--depth=2, chunkSize=100--depth=2.
