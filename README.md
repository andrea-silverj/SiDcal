# SiDcal
Similarity index calculator for codon usage analysis

usage: SiDcal.py [-h] [-a codon_table_A] [-b codon_table_B] [-ex] [-ref]

Hello dear User, and welcome to SiDcal!
This program calculates the similarity between two codon usage tables,
using the method described in the paper of Zhou et al., 2013

Interpretation of the values:
if codon_table_A=codon_table_B, SiD is equal to 0

optional arguments:
  -h, --help        show this help message and exit
  -a codon_table_A  |e.g., a virus codon table
  -b codon_table_B  |e.g., a host codon table
  -ex               |Print an example of the input format
  -ref              |Print the full reference of the Zhou et al.'s paper
