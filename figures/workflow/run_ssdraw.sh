#! /bin/bash

mkdir output
ssdraw single -f input/aligned.fasta -p input/2n1v.pdb -n 2n1v -o output/2n1v_white --fontcolor "red" --fontsize 12
ssdraw single -f input/aligned.fasta -p input/2n1v.pdb -n 2n1v -o output/2n1v_cmap -conservation_score --fontcolor "red" --fontsize 12
