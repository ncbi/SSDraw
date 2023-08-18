#! /bin/bash

mkdir output

python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_white

# conservation scores
python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_conservation_score_full-length -conservation_score
python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_conservation_score_60-225 -conservation_score --start 60 --end 225
python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_conservation_score_75-160 -conservation_score --start 75 --end 160

python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_bfactor_full-length -bfactor --color_map plasma
python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_bfactor_60-225 -bfactor --color_map plasma --start 60 --end 225
python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_bfactor_75-160 -bfactor --color_map plasma --start 75 --end 160

python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_mview_full-length -mview
python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_mview_60-225 -mview --start 60 --end 225
python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_mview_75-160 -mview --start 75 --end 160

python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_custom_full-length --scoring_file 1oqy_score_mod.txt --color_map viridis
python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_custom_60-225 --scoring_file 1oqy_score_mod.txt --color_map viridis --start 60 --end 225
python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_custom_75-160 --scoring_file 1oqy_score_mod.txt --color_map viridis --start 75 --end 160

python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_mview_nonsense -mview --color_map viridis
python3 ../../EAC_SSDraw.py -f aligned.fasta -p 1oqy.pdb -n 1oqy -o output/hHR23A_1oqy_solid_blue_nonsense --color_map viridis


