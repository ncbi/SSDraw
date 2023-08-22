#! /bin/bash

# example 1
python3 ../SSDraw.py --fasta 1ndd.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_out

# example 2: score by conservation
python3 ../SSDraw.py --fasta aligned.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_conservation -conservation_score

# example 3: score by bfactor
python3 ../SSDraw.py --fasta 1ndd.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_bfactor -bfactor

# example 4: Custom scoring file with custom color map
python3 ../SSDraw.py --fasta 2kdl.fasta --name 2kdl --pdb 2kdl.pdb --output 2kdl_out --scoring_file 2kdl_scoring.txt --color_map black cyan  

#Example 5: Choose subregion of alignment to run SSDraw on
python3 ../SSDraw.py --fasta aligned.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_conservation_cropped -conservation_score --start 80 --end 132
