#! /bin/bash

python3 ../../SSDraw.py --fasta 5ond.fasta --name 5OND_A --pdb 5ond.pdb --output 5ond_ungapped --color '#4C03E0'

python3 ../../SSDraw.py --fasta 5ond_gaps.fasta --name 5OND_A --pdb 5ond.pdb --output 5ond_gapped --color '#4C03E0'

python3 ../../SSDraw.py --fasta RfaH_alignment.fasta --name 5OND_A --pdb 5ond.pdb --output 5ond_aligned --color '#4C03E0'

python3 ../../SSDraw.py --fasta RfaH_alignment.fasta --name NusG_C --pdb NusG.pdb --output NusG_aligned --color '#9C6FF7' --chain_id C --SS NusG.horiz
