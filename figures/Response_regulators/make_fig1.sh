#! /bin/bash

python3 ../../SSDraw.py --fasta 5XSO_4KFC_alignment.fa --name 4KFC_1 --pdb 4kfc_bfac.pdb --output 4kfc_aligned --scoring_file 4kfc_score.txt --color_map gray red

python3 ../../SSDraw.py --fasta 4kfc.fasta --name 4KFC_1 --pdb 4kfc_bfac.pdb --output 4kfc_unaligned --scoring_file 4kfc_score.txt --color_map gray red

python3 ../../SSDraw.py --fasta 5XSO_4KFC_alignment.fa --name 5XSO_1 --pdb 5xso_bfac.pdb --output 5xso_aligned --scoring_file 5xso_score.txt --color_map gray blue

python3 ../../SSDraw.py --fasta 5xso.fasta --name 5XSO_1 --pdb 5xso_bfac.pdb --output 5xso_unaligned --scoring_file 5xso_score.txt --color_map gray blue


