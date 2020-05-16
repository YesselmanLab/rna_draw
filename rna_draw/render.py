import os
import time
import argparse
import cairosvg

from rna_draw import render_rna, svg, inv_utils, parameters
from rna_draw.colorer import *


import re
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-seq', help='seqs', required=False)
    parser.add_argument('-ss', help='ss', required=True)
    args = parser.parse_args()
    return args


def draw_rna(sequence, secstruct, colors, filename, parameters):
    r = render_rna.RNARenderer()

    pairmap = render_rna.get_pairmap_from_secstruct(secstruct)
    pairs = []
    for i in range(len(pairmap)):
        if pairmap[i] > i:
            pairs.append({"from":i, "to":pairmap[i], "p":1.0, "color": COLORS["e"]})
    r.setup_tree(secstruct, parameters.NODE_R, parameters.PRIMARY_SPACE, parameters.PAIR_SPACE)
    size = r.get_size()

    cell_size = max(size) + parameters.CELL_PADDING * 2

    svgobj = svg.svg("%s.svg" % filename, cell_size, cell_size)
    r.draw(svgobj, parameters.CELL_PADDING, parameters.CELL_PADDING, colors, pairs,
           sequence, parameters.RENDER_IN_LETTERS)


class RNARender(object):
    def __init__(self):
        pass

    def render(self, ss, seq=None, colors=None, filename='secstruct', png=True):
        final_color_rbgs = []

        if type(colors) is str:
            if len(colors) != len(ss):
                raise ValueError(
                        "color string was specified but was not the same length as ss")
            final_color_rbgs = parse_color_single_letter_codes(colors)

        elif type(colors) is list:
            final_color_rbgs = colors

        elif colors is None:
            final_color_rbgs = color_by_res_type(seq, ss)

        else:
            print(colors)
            raise ValueError("unknown color type")

        p = parameters.Parameters(seq, ss)
        draw_rna(seq, ss, final_color_rbgs, filename, p)

        if png:
            cairosvg.svg2png(url=filename+".svg", write_to=filename+".png")
        
        #if png:
         #   os.system(
          #      '/Applications/Inkscape.app/Contents/Resources/bin/inkscape --export-area-drawing --export-dpi 150 --export-png $(pwd)/%s.png $(pwd)/%s.svg > /dev/null' % (
           #         filename, filename))
           # os.remove(filename+".svg")




def main():
    rr = RNARender()
    rr.render("((((....))))","GGGGAAAACCCC")
    exit()
    args = parse_args()

    if not args.seq:
        args.seq = " " * len(args.ss)

    colors = color_by_res_type(args.seq)
    args.seq = "".join(args.seq)
    filename = "secstruct"
    p = parameters.Parameters()

    draw_rna(args.seq, args.ss, colors, filename, p)
    os.system(
        '/Applications/Inkscape.app/Contents/Resources/bin/inkscape --export-dpi 600 --export-png $(pwd)/%s.png $(pwd)/%s.svg' % (
        "secstruct", "secstruct"))
    os.remove("secstruct.svg")



if __name__ == "__main__":
    main()



















