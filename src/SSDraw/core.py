import re
import sys
import argparse
import warnings
import numpy as np
import typing as T
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.colors as mcolors

from Bio import pairwise2
from Bio import AlignIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.PDB.DSSP import DSSP
from Bio.Align import substitution_matrices
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from matplotlib.colors import ListedColormap


SPACING = 2.3


def read_r4s(input_file: str) -> T.Tuple[str, T.List[int]]:
    # algorithm to convert raw scores to grades is taken from consurf:
    # https://github.com/Rostlab/ConSurf
    seq = ""
    scores = []

    pattern = r"^\s*?\d+\s+(\w)\s+(\S+)\s+\[\s*\S+,\s*\S+\]\s+\S+\s+\d+\/\d+"
    with open(input_file) as READ:
        for line in READ:
            line = line.strip()
            if re.match(pattern, line):
                match = re.match(pattern, line)
                seq += match.group(1)
                scores.append(eval(match.group(2)))

    max_cons = min(scores)

    # 9 steps from -|max_cons| to |max_cons|, midpoint is 0
    ConsGradeUnity = max_cons / 4.5 * -1
    if max_cons >= 0:
        ConsGradeUnity = max_cons

    grades = []
    for score in scores:
        grades.append(max(1, 9 - int((score - max_cons) / ConsGradeUnity)))

    return seq, grades


def read_consurf_grad(input_file: str) -> T.Tuple[str, T.List[int]]:
    grades = []
    seq = ""
    pattern = r"^\s*?\d+\s+(\w)\s+\S+\s+\S+\s+(\d)\S*\s+-?\d+.\d+,\s+-?\d+.\d+\s+\d,\d\s+\d+\/\d+\s+\S+"
    with open(input_file, "r") as f:
        for line in f:
            if re.match(pattern, line):
                match = re.match(pattern, line)
                seq += match.group(1)
                grades.append(match.group(2))

    return seq, grades


def check_consurf_file(file: str) -> T.Optional[str]:
    consurf_pattern = r"^\s*?\d+\s+(\w)\s+\S+\s+\S+\s+(\d)\S*\s+-?\d+.\d+,\s+-?\d+.\d+\s+\d,\d\s+\d+\/\d+\s+\S+"
    r4s_pattern = (
        r"^\s*?\d+\s+(\w)\s+(\S+)\s+\[\s*\S+,\s*\S+\]\s+\S+\s+\d+\/\d+"
    )
    with open(file, "r") as f:
        for line in f:
            if re.match(consurf_pattern, line):
                return "consurf"
            if re.match(r4s_pattern, line):
                return "r4s"


def gap_sequence(seq: T.Any, extra_gaps: T.List[int]) -> T.Any:
    # seq can be a list or a string, anything that can be indexed
    # extra gaps is a list of length two [x,y],
    # where x is the number of characters to remove from the beginning
    # and y the number of characters to remove from the end
    new_seq = seq
    if extra_gaps[1] != 0:
        new_seq = new_seq[: -extra_gaps[1]]
    return new_seq[extra_gaps[0] :]


def NormalizeData(data: np.ndarray) -> np.ndarray:
    if np.min(data) == np.max(data):
        warnings.warn("Warning: scores are the same for all residues")
        return data
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def coords2path(
    coord_set1: T.List[T.Any],
) -> T.Tuple[T.List[T.Any], T.List[int]]:
    coords_f1 = []
    instructions1 = []

    for c in coord_set1:
        for n in range(len(c)):
            coords_f1.append(c[n])
            if n == 0:
                instructions1.append(1)
            else:
                instructions1.append(2)

    return coords_f1, instructions1


def build_loop(
    loop: T.Tuple[int, int],
    idx: int,
    ssidx: float,
    loop_coords: T.List[np.ndarray],
    linelen: int,
    nlines: int,
    prev_ss: T.Optional[str],
    next_ss: T.Optional[str],
    z: int = 1,
    clr: str = "r",
    mat: int = 0,
    size: int = 75,
) -> None:
    # if loop is smaller than 3 residues and has gaps on both sides, don't draw
    if prev_ss == "B" and next_ss == "B" and loop[1] - loop[0] < 2:
        return

    i0 = loop[0]
    if loop[0] != 0 and prev_ss != "B":
        i0 = loop[0] - 1
    elif i0 == 0:
        i0 = 0.06
    else:
        i0 = loop[0]
    i1 = loop[1] + 2
    if loop[1] == linelen - 1:
        i1 += 2

    o = 2
    if idx == nlines - 1:
        o = 0
    if next_ss == "B":
        o = -1.5
    if next_ss == None:
        o = -4.1

    rectangle = mpatch.Rectangle(
        (i0 / 6.0, -0.25 - 5.5 * idx - SPACING * ssidx),
        (i1 - i0 + o) / 6.0,
        0.5,
        fc=clr,
        ec="k",
        zorder=0,
    )

    xy = rectangle.get_xy()
    w = rectangle.get_width()
    h = rectangle.get_height()
    loop_coords.append(
        np.array(
            [
                [xy[0], xy[1]],
                [xy[0], xy[1] + h],
                [xy[0] + w, xy[1] + h],
                [xy[0] + w, xy[1]],
                [xy[0], xy[1]],
            ]
        )
    )


def build_strand(
    strand: T.Tuple[int, int],
    idx: int,
    ssidx: float,
    strand_coords: T.List[T.Any],
    next_ss: T.Optional[str],
    z: int = 1,
    clr: str = "r",
    imagemat: int = 0,
    size: int = 75,
) -> None:
    delta = 0 if next_ss == None else 1

    arrow = mpatch.FancyArrow(
        ((strand[0] + delta - 1) / 6.0),
        -5.5 * idx - SPACING * ssidx,
        (strand[1] - strand[0] + 1) / 6.0,
        0,
        width=1.0,
        fc=clr,
        linewidth=0.5,
        ec="k",
        zorder=z,
        head_width=2.0,
        length_includes_head=True,
        head_length=2.0 / 6.0,
    )

    strand_coords.append(arrow.get_xy())


def build_helix(
    helix: T.Tuple[int, int],
    idx: int,
    ssidx: float,
    coord_set1: T.List[T.Any],
    coord_set2: T.List[T.Any],
    clr: str = "r",
    size: float = 37.5,
    z: int = 1,
    bkg: T.Tuple[float, float, float] = (0.195, 0, 0.051),
    imagemat: int = 0,
) -> None:

    i = helix
    l = i[1] - i[0] + 1
    points = [
        [i[0] / 6.0, -0.25 - 5.5 * idx - SPACING * ssidx],
        [i[0] / 6.0 + 1.0 / 6, 0.75 - 5.5 * idx - SPACING * ssidx],
        [i[0] / 6.0 + 2.0 / 6, 0.75 - 5.5 * idx - SPACING * ssidx],
        [i[0] / 6.0 + 1.0 / 6, -0.25 - 5.5 * idx - SPACING * ssidx],
    ]
    # hlx = plt.Polygon(points,fc=clr,ec='k',zorder=1,linewidth=2)
    # coords= hlx.get_xy()
    # coord_set2.append(coords)
    coord_set2.append(points + [points[0]])

    for j in range((l - 2) - 1):
        if j % 2 == 0:
            points = [
                [
                    i[0] / 6.0 + (1.0 + j) / 6,
                    0.75 - 5.5 * idx - SPACING * ssidx,
                ],
                [
                    i[0] / 6.0 + (2.0 + j) / 6,
                    0.75 - 5.5 * idx - SPACING * ssidx,
                ],
                [
                    i[0] / 6.0 + (3.0 + j) / 6,
                    -0.75 - 5.5 * idx - SPACING * ssidx,
                ],
                [
                    i[0] / 6.0 + (2.0 + j) / 6,
                    -0.75 - 5.5 * idx - SPACING * ssidx,
                ],
            ]
            coord_set1.append(points + [points[0]])
            # hlx = mpatch.Polygon(points,fc=bkg,zorder=z)

        else:
            points = [
                [
                    i[0] / 6.0 + (1.0 + j) / 6,
                    -0.75 - 5.5 * idx - SPACING * ssidx,
                ],
                [
                    i[0] / 6.0 + (2.0 + j) / 6,
                    -0.75 - 5.5 * idx - SPACING * ssidx,
                ],
                [
                    i[0] / 6.0 + (3.0 + j) / 6,
                    0.75 - 5.5 * idx - SPACING * ssidx,
                ],
                [
                    i[0] / 6.0 + (2.0 + j) / 6,
                    0.75 - 5.5 * idx - SPACING * ssidx,
                ],
            ]
            coord_set2.append(points + [points[0]])
            # hlx = mpatch.Polygon(points,fc=clr,zorder=0)

    if (l - 2 - 1) % 2 == 1:

        points = [
            [i[1] / 6.0 - 1.0 / 6, -0.75 - 5.5 * idx - SPACING * ssidx],
            [i[1] / 6.0, -0.75 - 5.5 * idx - SPACING * ssidx],
            [i[1] / 6.0 + 1.0 / 6, 0.25 - 5.5 * idx - SPACING * ssidx],
            [i[1] / 6.0, 0.25 - 5.5 * idx - SPACING * ssidx],
        ]

        coord_set2.append(points + [points[0]])
        # hlx = mpatch.Polygon(points,fc=clr,zorder=0)

    else:
        points = [
            [i[1] / 6.0 - 1.0 / 6, 0.75 - 5.5 * idx - SPACING * ssidx],
            [i[1] / 6.0, 0.75 - 5.5 * idx - SPACING * ssidx],
            [i[1] / 6.0 + 1.0 / 6, -0.25 - 5.5 * idx - SPACING * ssidx],
            [i[1] / 6.0, -0.25 - 5.5 * idx - SPACING * ssidx],
        ]
        coord_set1.append(points + [points[0]])

        # hlx = plt.Polygon(points,fc=bkg,zorder=10)


def SS_breakdown(
    ss: str,
) -> T.Tuple[
    T.List[T.Tuple[int, int]],
    T.List[T.Tuple[int, int]],
    T.List[T.Tuple[int, int]],
    T.List[T.Tuple[int, int]],
    T.List[str],
    T.List[T.Tuple[int, int]],
]:
    i = 0
    curSS = ""
    jstart = -1
    jend = -1

    strand = []
    loop = []
    helix = []
    ssbreak = []

    ss_order = []
    ss_bounds = []

    last_ss = ""

    SS_equivalencies = {
        "H": ["H"],
        "-": ["-"],
        "S": [" ", "S", "C", "T", "G", "I", "P"],
        " ": [" ", "S", "C", "T", "G", "I", "P"],
        "C": [" ", "S", "C", "T", "G", "I", "P"],
        "T": [" ", "S", "C", "T", "G", "I", "P"],
        "G": [" ", "S", "C", "T", "G", "I", "P"],
        "I": [" ", "S", "C", "T", "G", "I", "P"],
        "P": [" ", "S", "C", "T", "G", "I", "P"],
        "E": ["E", "B"],
        "B": ["E", "B"],
    }

    cur_SSDict = {
        "H": "helix",
        "E": "strand",
        "B": "strand",
        "-": "break",
    }

    for i in range(len(ss)):
        if i == 0:
            curSS = SS_equivalencies[ss[i]]
            jstart = i
            if ss[i] in cur_SSDict.keys():
                last_ss = cur_SSDict[ss[i]]
            else:
                last_ss = "loop"
            continue

        if ss[i] in curSS:
            jend = i

        if ss[i] not in curSS or i == len(ss) - 1:
            if "E" in curSS and jend - jstart + 1 >= 3:
                strand.append((jstart, jend))
                ss_bounds.append((jstart, jend))
                ss_order.append("E")
                last_ss = "strand"
            elif "H" in curSS and jend - jstart + 1 >= 4:
                helix.append((jstart, jend))
                ss_bounds.append((jstart, jend))
                ss_order.append("H")
                last_ss = "helix"
            elif " " in curSS and last_ss != "loop":
                if jend < jstart:
                    jend = jstart
                loop.append((jstart, jend))
                ss_bounds.append((jstart, jend))
                ss_order.append("L")
                last_ss = "loop"
            elif "-" in curSS:
                if jend < jstart:
                    jend = jstart
                ssbreak.append((jstart, jend))
                ss_bounds.append((jstart, jend))
                ss_order.append("B")
                last_ss = "break"
            elif last_ss == "loop":
                if jend < jstart:
                    jend = jstart
                if len(loop) > 0:
                    jstart = loop[-1][0]
                    loop = loop[0:-1]
                    ss_bounds = ss_bounds[0:-1]
                    ss_order = ss_order[0:-1]
                loop.append((jstart, jend))
                ss_bounds.append((jstart, jend))
                ss_order.append("L")
                last_ss = "loop"
            else:
                if jend < jstart:
                    jend = jstart
                loop.append((jstart, jend))
                ss_bounds.append((jstart, jend))
                ss_order.append("L")
                last_ss = "loop"

            jstart = i
            curSS = SS_equivalencies[ss[i]]

    return strand, loop, helix, ssbreak, ss_order, ss_bounds


def updateSS(ss: str, seq: str, alignment: str, ref_align: str) -> str:
    ss_u = ""
    j = 0
    for i in range(len(alignment)):
        if alignment[i] == "-":
            if ref_align[i] == "-":
                ss_u += "-"
            else:
                ss_u += "C"
        else:
            ss_u += ss[j]
            j += 1

    return ss_u


def SS_align(
    alignment: T.List[str],
    ID: str,
    seq: str,
    ss: str,
    i_start: int,
    i_end: int,
) -> T.Tuple[str, str, T.List[int], int, int]:
    a_seq = ""
    seq_found = 0

    for i in alignment:

        if seq_found and i[0] == ">":
            break

        if i[0] == ">" and bool(re.search(ID.lower(), i.lower())):
            seq_found = 1
            continue

        if seq_found and i[0] != ">":
            a_seq += i

    if i_end != 0:
        i_end = len(a_seq) - i_end - 1

    a_seq = gap_sequence(a_seq, [i_start, i_end])

    a = pairwise2.align.localxs(seq, a_seq, -1, -0.5)

    # check if the dssp annotation has any extra residues not in the alignment
    if a[0][1] != a_seq:
        print("extra residues in pdb found\n")

    # check how many gap marks are at the end and beginning of the alignment
    # (a_seq) and compare to the amount found in a[0][1]
    a_seq_gaps = [0, 0]
    new_aln_gaps = [0, 0]
    for i in range(len(a_seq)):
        if a_seq[i] == "-":
            a_seq_gaps[0] += 1
        else:
            break
    for i in range(len(a_seq) - 1, -1, -1):
        if a_seq[i] == "-":
            a_seq_gaps[1] += 1
        else:
            break

    for i in range(len(a[0][1])):
        if a[0][1][i] == "-":
            new_aln_gaps[0] += 1
        else:
            break
    for i in range(len(a[0][1]) - 1, -1, -1):
        if a[0][1][i] == "-":
            new_aln_gaps[1] += 1
        else:
            break

    extra_gaps = [
        new_aln_gaps[0] - a_seq_gaps[0],
        new_aln_gaps[1] - a_seq_gaps[1],
    ]

    SS_updated = updateSS(ss, seq, a[0][0], a[0][1])

    SS_updated_new = gap_sequence(SS_updated, extra_gaps)
    a_new = gap_sequence(a[0][1], extra_gaps)

    return SS_updated_new, a_new, extra_gaps, i_start, i_end


def plot_coords(
    coords_all: T.List[T.List[T.Any]],
    mat: np.ndarray,
    sz: float,
    CMAP: T.Any,
    plot: T.Optional[T.Any] = None,
    ysz: float = 0.5,
) -> None:
    for i, coords in enumerate(coords_all):

        if not coords:
            continue

        coords_f1, instructions1 = coords2path(coords)

        # If loop or bottom helix layer, zorder = 0
        if i in [0, 1]:
            z = 0
        else:
            z = 10
        path = mpath.Path(np.array(coords_f1), np.array(instructions1))
        patch = mpatch.PathPatch(path, facecolor="none", ec="k", zorder=z)
        if plot != None:
            plot.add_patch(patch)
        else:
            plt.gca().add_patch(patch)
        im = plt.imshow(
            mat,
            extent=[0.0, sz, ysz, 3],
            cmap=CMAP,
            interpolation="none",
            zorder=z,
        )
        im.set_clip_path(patch)


def run_dssp(
    pdb_path: str, id: str, chain_id: str, dssp_exe: str = "mkdssp"
) -> T.List[str]:
    dssp_mode = "dssp"
    ss_seq = ""
    aa_seq = ""
    # Suppress Biopython PDBParser warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", PDBConstructionWarning)
        p = PDBParser()
        structure = p.get_structure(id, pdb_path)
    model = structure[0]

    try:
        dssp = DSSP(model, pdb_path, dssp=dssp_exe)
    except:
        # use pydssp instead of dssp
        dssp_mode = "pydssp"
        import torch
        import pydssp

        atom_indices = ["N", "CA", "C", "O"]
        # res_len = len([])

        for chain in model:
            if chain.get_id() == chain_id:
                resnum = len([i for i in chain.get_residues()])
                coords = torch.zeros(resnum, 4, 3)
                r_idx = 0
                for residue in chain:
                    for atom in residue:
                        if atom.name in atom_indices:
                            coords[r_idx][atom_indices.index(atom.name)] = (
                                torch.from_numpy(atom.coord)
                            )
                    aa_seq += seq1(residue.get_resname())
                    r_idx += 1

        ss_seq = "".join(pydssp.assign(coords, out_type="c3"))
        ss_seq = ss_seq.replace("-", "C")

    if dssp_mode == "dssp":
        a_key = list(dssp.keys())
        for key in a_key:
            if key[0] == chain_id:
                aa_seq += dssp[key][1]
                if dssp[key][2] == "-":
                    ss_seq += "C"
                else:
                    ss_seq += dssp[key][2]

    return [ss_seq, aa_seq]


def convert2horiz(dssp_file: str, pdbseq: str) -> T.List[str]:
    ss_seq = ""
    aa_seq = ""

    with open(dssp_file, "r") as f:
        lines = f.readlines()

    dssp_ext = dssp_file.split(".")[-1]

    if dssp_ext == "horiz":
        ss_seq = lines[0].rstrip("\n")
        aa_seq = lines[1].rstrip("\n")

    elif dssp_ext == "dssp":
        start_read = False
        for line in lines:
            if start_read:
                if line[13] == "!":
                    start_read = False
                    continue
                ss_seq += line[16]
                aa_seq += line[13]
            if line.split()[0] == "#":
                start_read = True

    else:
        raise Exception(
            "DSSP file extension not recognized: must be in .dssp or .horiz format"
        )

    return [ss_seq, aa_seq]


def score_column(msa_col: T.List[str], threshold: int = 0) -> float:
    blosum62 = substitution_matrices.load("BLOSUM62")
    # find consensus of the column
    aa_count = {
        "A": 0,
        "R": 0,
        "N": 0,
        "D": 0,
        "C": 0,
        "Q": 0,
        "E": 0,
        "G": 0,
        "H": 0,
        "I": 0,
        "L": 0,
        "K": 0,
        "M": 0,
        "F": 0,
        "P": 0,
        "S": 0,
        "T": 0,
        "W": 0,
        "Y": 0,
        "V": 0,
    }
    for i in msa_col:
        try:
            aa_count[i] += 1
        except:
            pass
    consensus_aa = max(zip(aa_count.values(), aa_count.keys()))[1]

    conservation_count = 0
    for i in msa_col:
        if i in aa_count.keys():
            if blosum62[consensus_aa][i] >= 0:
                conservation_count += 1

    return conservation_count / len(msa_col)


def parse_color(
    args: argparse.Namespace,
    seq_wgaps: str,
    pdbseq: str,
    bfactors: T.List[float],
    msa: T.List[str],
    extra_gaps: T.List[int],
) -> T.Tuple[T.Any, T.List[T.Any]]:
    CMAP = ""
    if (
        args.color in mcolors.BASE_COLORS.keys()
        or args.color in mcolors.CSS4_COLORS.keys()
        or args.color in mcolors.XKCD_COLORS.keys()
    ):
        CMAP = ListedColormap([args.color])
    elif args.color[0] == "#":
        CMAP = ListedColormap([args.color])
    if args.conservation_score or args.bfactor or args.scoring_file:
        if len(args.color_map) == 1:
            CMAP = args.color_map[0]
        else:
            CMAP = ListedColormap(args.color_map)

    # bvals are to make the colormap; taken from input PDB
    bvals = []

    if args.mview:
        mview_colors = {
            "A": 0,
            "G": 0,
            "I": 0,
            "L": 0,
            "M": 0,
            "P": 0,
            "V": 0,
            "F": 1,
            "H": 1,
            "W": 1,
            "Y": 1,
            "K": 2,
            "R": 2,
            "D": 3,
            "E": 3,
            "S": 4,
            "T": 4,
            "N": 5,
            "Q": 5,
            "C": 6,
        }
        mview_colors_hit = [0, 0, 0, 0, 0, 0, 0]

        mview_color_map = [
            "#33cc00",
            "#009900",
            "#cb0000",
            "#0133ff",
            "#0299fe",
            "#6601cc",
            "#ffff00",
            "#808080",
        ]

        for i in range(len(seq_wgaps)):
            try:
                m = mview_colors[seq_wgaps[i]]
                bvals.append(m)
                mview_colors_hit[m] += 1
            except:
                bvals.append(7)

        # remove colors of residues not in sequence
        for i in range(len(mview_colors_hit)):
            if mview_colors_hit[i] == 0:
                mview_color_map.pop(i)
                for j in range(len(bvals)):
                    if bvals[j] > i:
                        bvals[j] -= 1

        CMAP = ListedColormap(mview_color_map)

    elif args.consurf:
        # read in a rate4site output file or consurf score file
        # if rate4site output file, convert raw scores to 1-9 scores
        consurf_color_map = [
            "#10C8D1",
            "#8CFFFF",
            "#D7FFFF",
            "#EAFFFF",
            "#FFFFFF",
            "#FCEDF4",
            "#FAC9DE",
            "#F07DAB",
            "#A02560",
        ]

        scoring_seq = ""
        bvals_tmp = []
        consurf_mode = check_consurf_file(args.consurf)
        if consurf_mode == "consurf":
            scoring_seq, bvals_tmp = read_consurf_grad(args.consurf)
        elif consurf_mode == "r4s":
            scoring_seq, bvals_tmp = read_r4s(args.consurf)

        # print(bvals_tmp)

        # remove colors of residues not in sequence
        for i in reversed(range(1, 10)):
            if i not in bvals_tmp:
                consurf_color_map.pop(i - 1)
                for j in range(len(bvals_tmp)):
                    if bvals_tmp[j] > i:
                        bvals_tmp[j] -= 1
        # print(consurf_color_map)

        CMAP = ListedColormap(consurf_color_map)

        score_align = pairwise2.align.localxs(pdbseq, scoring_seq, -1, -0.5)

        j = 0
        for i in range(len(score_align[0][1])):
            if score_align[0][0][i] != "-":
                if score_align[0][1][i] != "-":
                    bvals.append(bvals_tmp[j])
                    j += 1
                else:
                    bvals.append(min(bvals_tmp))

    elif args.scoring_file:  # use custom scoring by residue
        # read in scoring file
        bvals_tmp = []
        scoring_seq = ""
        with open(args.scoring_file, "r") as g:
            lines = g.readlines()

        for line in lines:
            scoring_seq += line.split()[0]
            bvals_tmp.append(float(line.split()[1]))

        score_align = pairwise2.align.localxs(pdbseq, scoring_seq, -1, -0.5)

        j = 0
        for i in range(len(score_align[0][1])):
            if score_align[0][0][i] != "-":
                if score_align[0][1][i] != "-":
                    bvals.append(bvals_tmp[j])
                    j += 1
                else:
                    bvals.append(min(bvals_tmp))

    elif args.bfactor:  # score by bfactor
        bvals = [b for b in bfactors]

    elif args.conservation_score:  # score by conservation score
        bvals = []
        for i in range(len(msa[0])):
            bvals.append(score_column([msa[j][i] for j in range(len(msa))]))

    else:  # solid color
        bvals = [i for i in range(len(msa[0]))]

    if len(bvals) == len(pdbseq):
        # remove extra residues
        pdbseq = gap_sequence(pdbseq, extra_gaps)
        bvals = gap_sequence(bvals, extra_gaps)

        j = 0
        bvalsf = []
        bvals_align = pairwise2.align.localxs(seq_wgaps, pdbseq, -1, -0.5)
        if bvals_align[0][0] != seq_wgaps:
            print("Error in alignment or pdb sequence")
            sys.exit()
        for i in range(len(seq_wgaps)):
            if seq_wgaps[i] == "-" or bvals_align[0][1][i] == "-":
                if j >= len(bvals):
                    bvalsf.append(bvals[j - 1])
                else:
                    bvalsf.append(bvals[j])
            elif seq_wgaps[i] == bvals_align[0][1][i]:
                bvalsf.append(bvals[j])
                j += 1

        bvals = bvalsf

    return CMAP, bvals


def read_pdb(id: str, args: argparse.Namespace) -> T.Tuple[T.List[float], str]:
    pdbseq = ""
    # Suppress Biopython PDBParser warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", PDBConstructionWarning)
        p = PDBParser()
        bfactors = []
        structure = p.get_structure(id, args.pdb)
        model = structure[0]
        for chain in model:
            if chain.get_id() == args.chain_id:
                for residue in chain:
                    for atom in residue:
                        if atom.name == "CA":
                            bfactors.append(atom.bfactor)
                            pdbseq += seq1(residue.get_resname())
                            break
    sys.stdout = sys.__stdout__

    return bfactors, pdbseq


def initialize(
    args: T.Optional[argparse.Namespace] = None,
    parser: T.Optional[argparse.ArgumentParser] = None,
) -> T.Tuple[
    argparse.Namespace,
    str,
    T.List[float],
    T.List[str],
    str,
    str,
    T.List[int],
    int,
    int,
    T.List[T.Tuple[int, int]],
    T.List[T.Tuple[int, int]],
    T.List[T.Tuple[int, int]],
    T.List[T.Tuple[int, int]],
    T.List[str],
    T.List[T.Tuple[int, int]],
]:
    id = args.name

    chain_id = args.chain_id
    print("\nRunning for: " + id)

    # read in amino acid sequence from PDB
    bfactors, pdbseq = read_pdb(id, args)

    if args.SS:
        # get secondary structure from pre-existing DSSP annotation
        f = convert2horiz(args.SS, pdbseq)
    else:
        # run the dssp executable
        f = run_dssp(args.pdb, id, chain_id, dssp_exe=args.dssp_exe)

    nlines = 1
    salign = open(args.fasta).read().splitlines()

    if args.start > args.end:
        raise Exception("--start cannot be greater than --end")

    #####Align secondary structure to match input sequence alignment
    ss_wgaps, seq_wgaps, extra_gaps, i_start, i_end = SS_align(
        salign, args.name, f[1], f[0], args.start, args.end
    )

    # Break down secondary structure classifications in to continuous
    # chunks of helix, strand, and coil
    strand, loop, helix, ss_break, ss_order, ss_bounds = SS_breakdown(ss_wgaps)

    msa = [
        gap_sequence(a, [i_start, i_end])
        for a in AlignIO.read(open(args.fasta), "fasta")
    ]

    return (
        args,
        pdbseq,
        bfactors,
        msa,
        ss_wgaps,
        seq_wgaps,
        extra_gaps,
        i_start,
        i_end,
        strand,
        loop,
        helix,
        ss_break,
        ss_order,
        ss_bounds,
    )


def SSDraw(
    args: T.Optional[argparse.Namespace] = None,
    parser: T.Optional[argparse.ArgumentParser] = None,
) -> None:
    (
        args,
        pdbseq,
        bfactors,
        msa,
        ss_wgaps,
        seq_wgaps,
        extra_gaps,
        i_start,
        i_end,
        strand,
        loop,
        helix,
        ss_break,
        ss_order,
        ss_bounds,
    ) = initialize(args)
    nlines = 1

    # Parse color and scoring args
    CMAP, bvals = parse_color(
        args, seq_wgaps, pdbseq, bfactors, msa, extra_gaps
    )

    mat = np.tile(NormalizeData(bvals), (100, 1))

    # set figure parameters
    sz = 0
    c = "none"
    bc = "none"

    # set sizes of SS chunks
    ss_prev = 0
    for i in range(len(ss_order)):

        if ss_order[i] == "H":
            ss_prev = ss_bounds[i][1] / 6.0 + 1 / 6.0
        else:
            ss_prev = ss_bounds[i][1] / 6.0

    if ss_order[-1] == "H":
        sz = ss_bounds[-1][1] / 6.0 + 1 / 6.0
    elif ss_order[-1] in ["E", "B"]:
        sz = ss_bounds[-1][1] / 6.0
    elif ss_order[-1] == "L":
        sz = (ss_bounds[-1][1]) / 6.0

    # Plot secondary structure chunks
    strand_coords = []
    loop_coords = []
    helix_coords1 = []
    helix_coords2 = []

    fig, ax = plt.subplots(ncols=1, figsize=(25, 2 + 1.5 * (nlines - 1)))

    for i in range(len(ss_order)):
        prev_ss = None
        next_ss = None
        if i != 0:
            prev_ss = ss_order[i - 1]
        if i != len(ss_order) - 1:
            next_ss = ss_order[i + 1]

        if ss_order[i] == "L":
            build_loop(
                ss_bounds[i],
                0,
                -(2.0 / SPACING),
                loop_coords,
                len(ss_wgaps),
                1,
                prev_ss,
                next_ss,
                z=0,
                clr=c,
                mat=mat,
                size=sz,
            )
        elif ss_order[i] == "H":
            build_helix(
                ss_bounds[i],
                0,
                -(2.0 / SPACING),
                helix_coords1,
                helix_coords2,
                z=i,
                clr=c,
                bkg=bc,
                imagemat=mat,
                size=sz,
            )
        elif ss_order[i] == "E":
            build_strand(
                ss_bounds[i],
                0,
                -(2.0 / SPACING),
                strand_coords,
                next_ss,
                z=i,
                clr=c,
                imagemat=mat,
                size=sz,
            )

    plot_coords(
        [loop_coords, helix_coords2, strand_coords, helix_coords1],
        mat,
        sz,
        CMAP,
    )

    seq_to_show = seq_wgaps

    for i, aa in enumerate(seq_to_show):
        ax.text(
            i / 6.0,
            0.2,
            aa,
            ha="center",
            va="bottom",
            fontsize=args.fontsize,
            fontfamily="monospace",
            color=args.fontcolor,
        )

    plt.ylim([0.5, 3])

    # remove spines and yticks
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
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
            i += res_x * args.ticks
            label_i += args.ticks
        ax.get_xaxis().set_ticks(ticks, labels=labels)
        ax.xaxis.set_ticks_position("top")

    ax.set_aspect(0.5)

    print(
        "Saving output to {:}.{:}...".format(
            args.output, args.output_file_type
        )
    )
    plt.savefig(
        args.output + "." + args.output_file_type,
        bbox_inches="tight",
        dpi=args.dpi,
        transparent=True,
    )
