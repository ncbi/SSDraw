![header](imgs/1ndd_conservation.png)
# SSDraw
SSDraw is a program that generates publication-quality protein secondary structure diagrams from three-dimensional protein structures. To depict relationships between secondary structure and other protein features, diagrams can be colored by conservation score, B-factor, or custom scoring.

SSDraw also has a colab notebook available at https://colab.research.google.com/github/ethanchen1301/SSDraw/blob/main/SSDraw.ipynb (only usable for Chrome)
## Installation

SSDraw requires the Biopython and matplotlib modules to be installed.

```bash
pip install biopython
pip install matplotlib
pip install numpy
```

SSDraw also requires DSSP to be installed in order to generate secondary structure annotations:

```bash
conda install -c salilab dssp
```

Alternatively, you can install DSSP either through apt-get (sudo apt-get install dssp), or you can follow the instructions on their github page to make a local installation: 
https://github.com/cmbi/dssp.

If DSSP fails to install, there is a PyTorch implementation of DSSP that you can install as a workaround:

#### Alternative to installing DSSP

```
pip install torch
pip install pydssp
pip install nomkl
```

## Instructions
SSDraw requires 4 arguments:
1. --fasta: the file name sequence or alignment file in fasta format.
2. --name: the id of the sequence in the fasta file corresponding to your protein of interest.
3. --pdb: the file name of the pdb file of your protein
4. --output: the output file name to use


#### Example 1:
```
python3 ../SSDraw.py --fasta 1ndd.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_out
```
![Example 1](imgs/1ndd_out.png)

### Coloring options:
SSDraw uses a gradient to color each position in the alignment by a certain score. The user can choose which scoring system to use, and they can also choose which colormap.

### Scoring: 
-conservation_score: score each position in the alignment by conservation score.

-bfactor: score each residue in the pdb by B-factor.

-scoring_file: score each residue by a custom scoring file prepared by the user.

-mview: color each residue by the mview coloring system.

#### Example 2: score by conservation
```
python3 ../SSDraw.py --fasta aligned.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_conservation -conservation_score --start 80 --end 132
```
![Example 2](imgs/1ndd_conservation.png)
Note: for more on how the --start and --end options work, see [choosing a subregion](#choosing-a-subregion).

#### Example 3: Score by B-factor
```
python3 ../SSDraw.py --fasta 1ndd.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_bfactor -bfactor
```
![Example 3](imgs/1ndd_bfactor.png)
### Choosing a colormap:
The default colormap for SSDraw is inferno. The user can select one of the matplotlib library color maps or simply list a space delimited set of colors they'd like to use with the --color_map option. Alternatively, the user can select a single color with the --color option and SSDraw will use that color on the whole image.

#### Example 4: Custom scoring file with custom color map
```
python3 ../SSDraw.py --fasta 2kdl.fasta --name 2kdl --pdb 2kdl.pdb --output 2kdl_out --scoring_file 2kdl_scoring.txt --color_map black cyan  
```
![Example 4](imgs/2kdl_out.png)
### DSSP files:
Normally, SSDraw will generate a DSSP annotation from the PDB file, but if you have a DSSP file you would like to use, you can upload it and input the file name in Options.

### Choosing a subregion:
If you want SSDraw to draw only a portion of your alignment, you can specify the start and/or end points using the --start and --end options respectively. The argument for these options correspond to the index of the alignment position, not to the residue position numbers. See [example 2](#example-2-score-by-conservation).

### Running on multiple pdbs:
In order to rapidly generate multiple images with SSDraw, we recommend writing shell scripts comprised of commands like those shown in the above examples. For examples of such a shell script, see [here](SSDraw/figures/ubiquitin/run_ssdraw.sh).

### Helper Scripts
We now provide two helper scripts to assist the user in stacking multiple diagrams.
```
python3 run_multiple_pdbs_on_one_msa.py --input [input script] --output [output directory]
```
run_multiple_pdbs_on_one_msa.py will run SSDraw for multiple pdbs from a single multiple sequence alignment, saving the diagrams to a specified output directory. Finally, the script will create a composite stacked image of the diagrams. An example input script is shown in SSDraw/example_run.txt.

```
python3 combine_images.py [list of images]
```
This script creates a composite stacked image in png format of a list of SSDraw diagrams. To use this script, you will already have had to run SSDraw to generate individual diagrams. 


