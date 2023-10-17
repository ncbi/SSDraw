"""
Reads in helper script to run SSDraw for multiple PDBs and 
one multiple sequence alignment, then combines images into a composite image

To run, run the command 
"python run_multiple_pdbs_on_one_msa.py --input [input script] --output [output name]"
An example input script is shown in "example_run.txt"
"""
import sys
import os
import re
from PIL import Image
import argparse
from combine_images import combine_images


def get_args():
    parser_description="A helper script to run SSDraw for multiple PDBs from a single multiple sequence alignment."
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = parser_description,
        epilog="")
    parser.add_argument("-i", "--input", help="name of input script")
    parser.add_argument("-o","--output", help="name of output image")
    
    args = parser.parse_args()
    return args, parser

def main():
    
    args,parser = get_args()
    output_file_type = "png"

    ssdraw_params = {
        "FASTA": [],
        "PDB": [],
        "NAME": [],
        "OUTPUT": [],
        "ADDITIONAL_PARAMS": []
    }

    with open(args.input, "r") as f:
        lines = f.readlines()

    current_param = ""
    read_state = False

    for line in lines:
        words = line.split()

        if len(words) > 0:

            if words[0] in ssdraw_params.keys():
                current_param = words[0]
                continue

            if words[0] == "{":
                read_state = True
                continue
            
            if words[0] == "}":
                read_state = False
                current_param = ""
            
            if words[0][0] == "#":
                continue
            
            if bool(re.search("--output_file_type", line)):
                output_file_type = line.strip()[19:]

            if current_param != "" and read_state:
                ssdraw_params[current_param].append(line.strip())


            

    # check if pdbs, names, and outputs are the same length

    if len(ssdraw_params["PDB"]) != len(ssdraw_params["NAME"]) or len(ssdraw_params["PDB"]) != len(ssdraw_params["OUTPUT"]):
        raise Exception("Number of options in PDB, NAME, and OUTPUT sections must be the same")
    

    # make dir
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    fasta = ssdraw_params["FASTA"][0]
    additional_params = " ".join(ssdraw_params["ADDITIONAL_PARAMS"])
    imgs = []
    for i in range(len(ssdraw_params["PDB"])):

        pdb = ssdraw_params["PDB"][i]
        name = ssdraw_params["NAME"][i]
        output = args.output+"/"+ssdraw_params["OUTPUT"][i]

        ssdraw_command = "python3 SSDraw.py -f {:} -p {:} -n {:} -o {:} {:}".format(fasta, pdb, name, output, additional_params)
        os.system(ssdraw_command)
        imgs.append(Image.open(output+"."+output_file_type))
        
    print("Creating composite image {:}/{:}.{:}".format(args.output,args.output, output_file_type))
    combine_images(imgs).save("{:}/{:}.{:}".format(args.output,args.output,output_file_type))
    

if __name__ == '__main__':
    main()