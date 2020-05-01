import render_rna
import os
import time
import svg
import inv_utils
import argparse
import re
import numpy as np

NODE_R = 10
PRIMARY_SPACE = 20
PAIR_SPACE = 23

CELL_PADDING = 40
TEXT_SIZE = 50

RENDER_IN_LETTERS = False

COLORS = {#"r": [255, 0, 0],
          "r": [255, 102, 102],
          "g": [113, 188, 120],
          "b": [51, 153, 255],
          #"b": [0, 153, 255],
          "k": [1, 0, 0],
          "y": [255, 211, 0],
          "c": [0, 255, 255],
          "m": [255, 0, 255],
          "w": [255, 255, 255],
          "e": [150, 150, 150],
          "o": [231, 115, 0],
          "i": [51, 204, 204],
          "h": [51, 153, 255]}
          #"i": [0, 204, 153],
          #"h": [46, 184, 46]}

def draw_rna(sequence, secstruct, colors, filename="secstruct"):
    r = render_rna.RNARenderer()

    pairmap = render_rna.get_pairmap_from_secstruct(secstruct)
    pairs = []
    for i in range(len(pairmap)):
        if pairmap[i] > i:
            pairs.append({"from":i, "to":pairmap[i], "p":1.0, "color":COLORS["e"]})
    r.setup_tree(secstruct, NODE_R, PRIMARY_SPACE, PAIR_SPACE)
    size = r.get_size()

    cell_size = max(size) + CELL_PADDING * 2
    colors = [COLORS[x] for x in list(colors)]

    svgobj = svg.svg("%s.svg" % filename, cell_size, cell_size)
    r.draw(svgobj, CELL_PADDING, CELL_PADDING, colors, pairs, sequence, RENDER_IN_LETTERS)


def parse_colors(color_string):
    colorings = color_string.strip().split(",")
    colors = []
    for coloring in colorings:
        if "x" in coloring:
            [n, color] = coloring.split("x")
            colors += int(n) * [color]
        else:
            colors += [coloring]
    return colors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-seq")
    parser.add_argument("-ss")

    args = parser.parse_args()

    if not args.seq:
        args.seq = " " * len(args.ss)

    final_seq = []
    colors = []
    for e in args.seq:
        if e == "A":
            colors.append("y")
            final_seq.append("A")
        elif e == "C":
            colors.append("g")
            final_seq.append("C")
        elif e == "G":
            colors.append("r")
            final_seq.append("G")
        elif e == "U" or e == "T":
            colors.append("b")
            final_seq.append("U")
        else:
            colors.append("e")
            final_seq.append(" ")
    args.seq = "".join(final_seq)

    draw_rna(args.seq, args.ss, colors)
    os.system(
        '/Applications/Inkscape.app/Contents/Resources/bin/inkscape --export-dpi 600 --export-png $(pwd)/%s.png $(pwd)/%s.svg' % (
        "secstruct", "secstruct"))



if __name__ == "__main__":
    main()



















