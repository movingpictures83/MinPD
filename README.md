# MinPD
# Language: C++
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 2.0, GCC 4.8.4


Distance-based phylogenetic analysis and recombination detection (Buendia and Narasimhan, 2004).

Program takes as input a TXT file of parameters, one per line:
[align=0|1]
        1= MinPD aligns sequences pairwise (no global alignment)
        0= no alignment of sequences
[fragments=1|4|8]
        1= no recombination detection
        4= 4 fragment recombination detection
        8= 8 fragment recombination detection

It will then output two files:
[prefix].d.txt: Ancestor-descendant relationships
[prefix].out.txt: NJ trees, to build the evolutionary framework
