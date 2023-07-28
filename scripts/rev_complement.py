#!/usr/bin/env python
import sys

# Complementing DNA
# defining the DNA dictionary
dna = sys.argv[1]
compl = {
    "A" : "T",
    "T" : "A",
    "C" : "G",
    "G" : "C"
}

rev = ""
for i in dna:
   rev = rev + compl[i]

print(rev)