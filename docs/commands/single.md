# SSDraw single Command

This section documents the `ssdraw single` command.

```bash
usage: generate_docs.py single [-h] [-f FASTA] [-p PDB] [-n NAME] [-o OUTPUT] [--SS SS] [--chain_id CHAIN_ID] [--color_map [COLOR_MAP ...]]
                               [--scoring_file SCORING_FILE] [--color COLOR] [-conservation_score] [--output_file_type OUTPUT_FILE_TYPE] [-bfactor]
                               [-mview] [--dpi DPI] [--ticks TICKS] [--start START] [--end END] [--dssp_exe DSSP_EXE] [--consurf CONSURF]
                               [--fontsize FONTSIZE] [--fontcolor FONTCOLOR]

options:
  -h, --help            show this help message and exit
  -f, --fasta FASTA     (required) sequence/alignment file in fasta format
  -p, --pdb PDB         (required) pdb file
  -n, --name NAME       (required) id of the protein in the alignment file
  -o, --output OUTPUT   (required) name for output file
  --SS SS               secondary structure annotation in DSSP or .horiz format. If this option is not provided, SSDraw will compute secondary structure
                        from the given PDB file with DSSP.
  --chain_id CHAIN_ID   chain id to use in pdb. Defaults to chain A.
  --color_map [COLOR_MAP ...]
                        color map to use for heat map
  --scoring_file SCORING_FILE
                        custom scoring file for alignment
  --color COLOR         color for the image. Can be a color name (eg. white, black, green), or a hex code
  -conservation_score   score alignment by conservation score
  --output_file_type OUTPUT_FILE_TYPE
                        output file type. Options: png, ps, eps, tif, svg
  -bfactor              score by B-factor
  -mview                color by mview color map
  --dpi DPI             dpi to use for final plot
  --ticks TICKS         set ticks at every nth position
  --start START
  --end END
  --dssp_exe DSSP_EXE   The path to your dssp executable. Default: mkdssp
  --consurf CONSURF     consurf or rate4site file to color image with. If rate4site file is given, SSDraw will convert raw scores to grades.
  --fontsize FONTSIZE   font size for residue numbers
  --fontcolor FONTCOLOR
                        font color for residue numbers
```
