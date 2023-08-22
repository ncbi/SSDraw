#! /bin/bash

mkdir output
python3 ../../SSDraw.py -f input/aligned.fasta -p input/2n1v.pdb -n 2n1v -o output/2n1v_white 
python3 ../../SSDraw.py -f input/aligned.fasta -p input/2n1v.pdb -n 2n1v -o output/2n1v_cmap -conservation_score
