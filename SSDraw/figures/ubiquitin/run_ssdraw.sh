#! /bin/bash

mkdir output
python3 ../../SSDraw.py -f input/aligned.fasta -p input/1ubq.pdb -n 1ubq -o output/UBQ_1ubq -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/1ndd.pdb -n 1ndd -o output/NEDD8_1ndd -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/2n1v.pdb -n 2n1v -o output/SUMO_2n1v -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/2l7r.pdb -n 2l7r -o output/FUBI_2l7r -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/1z2m.pdb -n 1z2m -o output/UCRP_1z2m -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/1oqy.pdb -n 1oqy -o output/hHR23A_1oqy -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/1p1a.pdb -n 1p1a -o output/hHR23B_1p1a -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/2bwf.pdb -n 2bwf -o output/Dsk2_2bwf -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/1yqb.pdb -n 1yqb -o output/Ubiquilin_1yqb -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/4k95.pdb -n 4k95 -o output/Parkin_4k95 -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/1vcb.pdb -n 1vcb -o output/Elongin_1vcb -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/1wh3.pdb -n 1wh3 -o output/OASL_1wh3 -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/2dzi.pdb -n 2dzi -o output/GDX_2dzi -conservation_score --start 70 --end 165
python3 ../../SSDraw.py -f input/aligned.fasta -p input/1wx9.pdb -n 1wx9 -o output/BAT3_1wx9 -conservation_score --start 70 --end 165