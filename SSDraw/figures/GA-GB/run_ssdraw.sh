#! /bin/bash

mkdir output
python3 ../../SSDraw.py -f input/rcsb_pdb_2KDL.fasta -p input/2kdl.pdb -n 2kdl -o output/2kdl --scoring_file input/scoring_files/2kdl_scoring.txt --color_map black cyan 
python3 ../../SSDraw.py -f input/rcsb_pdb_2KDM.fasta -p input/2kdm.pdb -n 2kdm -o output/2kdm --scoring_file input/scoring_files/2kdm_scoring.txt --color_map black black 
python3 ../../SSDraw.py -f input/rcsb_pdb_2LHC.fasta -p input/2lhc.pdb -n 2lhc -o output/2lhc --scoring_file input/scoring_files/2lhc_scoring.txt --color_map black yellow  
python3 ../../SSDraw.py -f input/rcsb_pdb_2LHD.fasta -p input/2lhd.pdb -n 2lhd -o output/2lhd --scoring_file input/scoring_files/2lhd_scoring.txt --color_map black magenta  
python3 ../../SSDraw.py -f input/rcsb_pdb_2LHG.fasta -p input/2lhg.pdb -n 2lhg -o output/2lhg --scoring_file input/scoring_files/2lhg_scoring.txt --color_map black green 
python3 ../../SSDraw.py -f input/rcsb_pdb_2LHE.fasta -p input/2lhe.pdb -n 2lhe -o output/2lhe --scoring_file input/scoring_files/2lhe_scoring.txt --color_map black white white  
