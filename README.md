![header](imgs/1ndd_conservation.png)
# SSDraw
SSDraw is a program that generates publication-quality protein secondary structure diagrams from three-dimensional protein structures. To depict relationships between secondary structure and other protein features, diagrams can be colored by conservation score, B-factor, or custom scoring.

SSDraw also has a colab notebook available at https://colab.research.google.com/github/ethanchen1301/SSDraw/blob/main/SSDraw.ipynb (only usable for Chrome)

## Installation

The easiest way to install SSDraw is via pip:

```bash
pip install SSDraw
```

You can find the package on PyPi [here](https://pypi.org/project/SSDraw/)

## Instructions
SSDraw provides two subcommands:
1. single: Run SSDraw on a single FASTA/PDB pair
2. multi: Run SSDraw for multiple PDBs from one MSA

SSDraw single requires 4 arguments:
1. --fasta: the file name sequence or alignment file in fasta format.
2. --name: the id of the sequence in the fasta file corresponding to your protein of interest.
3. --pdb: the file name of the pdb file of your protein
4. --output: the output file name to use


#### Example 1:
```
ssdraw single --fasta 1ndd.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_out
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
ssdraw single --fasta aligned.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_conservation -conservation_score --start 80 --end 132
```
![Example 2](imgs/1ndd_conservation.png)
Note: for more on how the --start and --end options work, see [choosing a subregion](#choosing-a-subregion).

#### Example 3: Score by B-factor
```
ssdraw single --fasta 1ndd.fasta --name 1ndd --pdb 1ndd.pdb --output 1ndd_bfactor -bfactor
```
![Example 3](imgs/1ndd_bfactor.png)
### Choosing a colormap:
The default colormap for SSDraw is inferno. The user can select one of the matplotlib library color maps or simply list a space delimited set of colors they'd like to use with the --color_map option. Alternatively, the user can select a single color with the --color option and SSDraw will use that color on the whole image.

#### Example 4: Custom scoring file with custom color map
```
ssdraw single --fasta 2kdl.fasta --name 2kdl --pdb 2kdl.pdb --output 2kdl_out --scoring_file 2kdl_scoring.txt --color_map black cyan  
```
![Example 4](imgs/2kdl_out.png)
### DSSP files:
Normally, SSDraw will generate a DSSP annotation from the PDB file, but if you have a DSSP file you would like to use, you can upload it and input the file name in Options.

### Choosing a subregion:
If you want SSDraw to draw only a portion of your alignment, you can specify the start and/or end points using the --start and --end options respectively. The argument for these options correspond to the index of the alignment position, not to the residue position numbers. See [example 2](#example-2-score-by-conservation).

### ConSurf
We now provide the option to color by ConSurf grade using --consurf option:
```
ssdraw single --fasta 2n1v.fasta --name 2n1v --pdb 2n1v.pdb --output 2n1v_consurf --consurf 2n1v_nscores.txt --ticks 5
```
![ConSurf](imgs/2n1v_consurf.png)
The code can read in either the ConSurf grades file or a Rate4Site file with raw scores. The raw scores will be converted into grades according to the algorithm used by ConSurf: https://github.com/Rostlab/ConSurf

### Stacking diagrams
We now provide a `multi` option to assist the user in stacking multiple diagrams.
```
ssdraw multi --input [input script] --output [output directory]
```

This will run SSDraw for multiple pdbs from a single multiple sequence alignment, saving the diagrams to a specified output directory. Finally, this will create a composite stacked image of the diagrams. An example input script is shown in SSDraw/examples/example_run.txt.

## Running with Docker

You can use Docker to run SSDraw without installing dependencies on your system. However, you must have Docker Desktop for your OS installed. You can find that [here](https://www.docker.com/get-started/).

Pull the image from the Docker Registry [here](https://hub.docker.com/r/prameshsharma25/ssdraw):
```
docker pull prameshsharma25/ssdraw
```

### 1. Run SSDraw with your input files

Make sure your input files (FASTA, PDB, etc.) are in your current directory.  
Then run:

```sh
docker run --rm -v "$PWD":/app prameshsharma25/ssdraw single --fasta your.fasta --name your_id --pdb your.pdb --output output_name
```

- Replace `your.fasta`, `your_id`, `your.pdb`, and `output_name` with your actual file names and sequence ID.
- The output image will be saved in your current directory.

#### Example

```sh
docker run --rm -v "$PWD":/app prameshsharma25/ssdraw single --fasta examples/1ndd.fasta --name 1ndd --pdb examples/1ndd.pdb --output 1ndd_out --dpi 1000
```

**Note:**  
The `-v "$PWD":/app` option mounts your current directory into the container, so SSDraw can access your files and save outputs locally.  
You can pass any SSDraw command-line options as usual after the image

### 3. Run SSDraw with multiple PDBs
Following the format in examples/example_run.txt, you're able to utilize SSDraw across several PDBs which builds a composite image of secondary structures.

```sh
docker run --rm -v "$PWD":/app prameshsharma25/ssdraw multi -i run.txt -o output_name
```

#### Example

```sh
docker run --rm -v "$PWD":/app prameshsharma25/ssdraw multi -i examples/example_run.txt -o ubiquitin_stacked
```


