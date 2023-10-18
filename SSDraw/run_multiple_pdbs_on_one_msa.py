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
import SSDraw
from combine_images import combine_images
import matplotlib.pyplot as plt
import numpy as np

def SSDraw_layer(fastas,pdbs,names,output_names,output_dir,additional_params):

    nlines = len(pdbs)

    #Plot secondary structure chunks
    strand_coords = []
    loop_coords = []
    helix_coords1 = []
    helix_coords2 = []

    #Parameters for scalable figure dimensions
    minsize = 999
    maxsize = -999
    loop_ys = []

    for i in range(len(pdbs)):
        additional_params += f' -f {fastas[0]} -p {pdbs[i]} -n {names[i]} -o {output_names[i]}'
        args,parser = SSDraw.get_args(args=additional_params.split())
        args,pdbseq,bfactors,msa,ss_wgaps,seq_wgaps,extra_gaps,i_start,i_end,strand,loop,helix,ss_break,ss_order,ss_bounds=SSDraw.initialize(args,parser)

        #Parse color and scoring args
        CMAP, bvals = SSDraw.parse_color(args,seq_wgaps,pdbseq,bfactors,msa,extra_gaps)

        mat = np.tile(SSDraw.NormalizeData(bvals), (100,1))

        #set figure parameters
        sz = 0
        c = 'none'
        bc = 'none'

        #set sizes of SS chunks
        ss_prev = 0
        for j in range(len(ss_order)):
        
            if ss_order[j] == 'H':
                ss_prev = ss_bounds[j][1]/6.0+1/6.0
            else:
                ss_prev = ss_bounds[j][1]/6.0

        if ss_order[-1] == 'H':
            sz = ss_bounds[-1][1]/6.0+1/6.0
        elif ss_order[-1] in ['E','B']:
            sz = ss_bounds[-1][1]/6.0
        elif ss_order[-1] == 'L':
            sz = (ss_bounds[-1][1])/6.0


        for j in range(len(ss_order)):
            prev_ss = None
            next_ss = None
            if j != 0:
                prev_ss = ss_order[j-1]
            if j != len(ss_order)-1:
                next_ss = ss_order[j+1]

            if ss_order[j] == 'L':
                SSDraw.build_loop(ss_bounds[j],1,i,loop_coords,len(ss_wgaps),nlines,prev_ss,next_ss,z=0,clr=c,mat=mat,size=sz)
            elif ss_order[j] == 'H':
                SSDraw.build_helix(ss_bounds[j],1,i,helix_coords1,helix_coords2,z=i,clr=c,bkg=bc,imagemat=mat,size=sz)
            elif ss_order[j] == 'E':
                SSDraw.build_strand(ss_bounds[j],1,i,strand_coords,next_ss,z=i,clr=c,imagemat=mat,size=sz)

        loop_ys.append(loop_coords[-1][0][1])
        

    for j in range(len(loop_coords)):

        if loop_coords[j][0][1] < minsize:
            minsize = loop_coords[j][0][1]

        if loop_coords[j][1][1] > maxsize:
            maxsize = loop_coords[j][1][1]

    #fig, ax = plt.subplots(ncols=1, figsize=(sz*.7,(2+1.5*(nlines-1)*.7)))
    
    fig, ax = plt.subplots(ncols=1, figsize=(sz*.7,((maxsize-minsize)*.37)))
    maxlen = max([len(x) for x in output_names])
    
    for i in range(len(output_names)):
        ax.annotate(output_names[i],xy=(0,loop_ys[i]),xytext=(-0.2,loop_ys[i]),fontsize=14,ha='right')
    

    SSDraw.plot_coords([loop_coords,helix_coords2,strand_coords,helix_coords1],mat,sz,CMAP,plot=plt.gca(),ysz=minsize-0.75)

    plt.ylim([minsize-0.75,maxsize+0.75])

    # remove spines and yticks
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_ticks([])

    if args.ticks == 0:
        ax.get_xaxis().set_ticks([])
    else:
        res_x = 0.1646
        ticks = []
        labels = []
        i = 0
        label_i = 1
        while label_i <= len(bvals):
            ticks.append(i)
            labels.append(str(label_i))
            i+=res_x*args.ticks
            label_i+=args.ticks
        ax.get_xaxis().set_ticks(ticks,labels=labels)
        ax.xaxis.set_ticks_position('top')


    ax.set_aspect(0.5)
    #plt.tight_layout()


def get_args():
    parser_description="A helper script to run SSDraw for multiple PDBs from a single multiple sequence alignment."
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = parser_description,
        epilog="")
    parser.add_argument("-i", "--input", help="name of input script")
    parser.add_argument("-o","--output", help="name of output directory")
    
    args = parser.parse_args()
    return args, parser

def parse_params(args):

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

    fasta = ssdraw_params["FASTA"][0]
    additional_params = " ".join(ssdraw_params["ADDITIONAL_PARAMS"])

    return ssdraw_params, additional_params,output_file_type
    

def main():

    if len(sys.argv) < 3:
        sys.exit("""usage: run_multiple_pdbs_on_one_msa.py [-h] [-i INPUT] [-o OUTPUT]

        A helper script to run SSDraw for multiple PDBs from a single multiple sequence alignment.

        options:
             -h, --help            show this help message and exit
             -i INPUT, --input INPUT
                        name of input script
             -o OUTPUT, --output OUTPUT
                        name of output directory""")
    
    args,parser = get_args()
    
    ssdraw_params,additional_params,output_file_type = parse_params(args)

    SSDraw_layer(ssdraw_params["FASTA"],ssdraw_params["PDB"],ssdraw_params["NAME"],ssdraw_params["OUTPUT"],args.output,additional_params)

    print("Creating composite image {:}.{:}".format(args.output, output_file_type))
    plt.savefig("{:}.{:}".format(args.output,output_file_type))

if __name__ == '__main__':
    main()
