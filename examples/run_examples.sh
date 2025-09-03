#! /bin/bash

# example 1
python3 ../SSDraw.py --fasta 1ndd.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_out --dpi 1000

# example 2: score by conservation
python3 ../SSDraw.py --fasta aligned.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_conservation -conservation_score --start 80 --end 162 --dpi 1000

# example 3: score by bfactor
python3 ../SSDraw.py --fasta 1ndd.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_bfactor -bfactor --dpi 1000

# example 4: Custom scoring file with custom color map
python3 ../SSDraw.py --fasta 2kdl.fasta --name 2kdl --pdb 2kdl.pdb --output 2kdl_out --scoring_file 2kdl_scoring.txt --color_map black cyan --dpi 1000

#example 5: automatically generate stacked diagrams
python3 ../run_multiple_pdbs_on_one_msa.py -i example_run.txt -o ubiquitin_stacked

#example 6: color by Rate4site scores and add ticks
python3 ../SSDraw.py -p 2n1v.pdb -f 2n1v.fasta -n 2n1v -o 2n1v_consurf --consurf 2n1v_nscores.txt --ticks 5
