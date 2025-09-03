#! /usr/local/bin/python3

import sys

aa3to1 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "MSE": "X",
    "NAG": "",
    "BMA": "",
    "MAN": "",
    "SO4": "",
    "NHE": "",
    " CL": "",
    "HOH": "",
    " CU": "",
    " ZN": "",
    "BNG": "",
    "PLM": "",
    "ACT": "",
    "ATM": "",
    " MG": "",
    "UNK": "X",
}

if __name__ == "__main__":

    f = open(sys.argv[1])
    startone = int(sys.argv[2])

    for i, l in enumerate(f):

        if l[:4] == "ATOM":
            if i > startone:
                print(aa3to1[l[17:20]], 1)
            else:
                print(aa3to1[l[17:20]], 0)
