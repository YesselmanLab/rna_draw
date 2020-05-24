import os
import time
import argparse
import cairosvg
import pandas as pd

from rna_draw import render_rna, svg, parameters, colorer
from rna_draw.colorer import *

import matplotlib.cm

import re
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-ss', help='secondary structure in dot bracket notation',
                        required=True)
    parser.add_argument('-seq', help='rna sequence', required=False)
    parser.add_argument('-out', help='output png file',
                        required=False, default='secstruct')
    parser.add_argument('-color_str', help='description of coloring, see docs for options',
                        required=False)
    parser.add_argument('-data_file', help='path to data to color by', required=False)
    parser.add_argument('-data', help='data values by res seperated by ;',required=False)
    parser.add_argument('-render_type', help='scheme to color by: res_type,paired,motif,none',
                        required=False, default="res_type")
    args = parser.parse_args()
    return args


class RNADrawer(object):
    def __init__(self):
        self.__colorer = colorer.Colorer()
        self.__draw_params = parameters.DrawParameters()


    def draw(self, ss, seq=None, filename='secstruct', color_str=None,
             render_type=None, data=None, data_palette=None, default_color=None,
             draw_params=None):

        self.__setup(draw_params)

        if seq is None:
            seq = " " * len(ss)

        final_color_rbgs = self.__get_rgb_colors(
                ss, seq, color_str, data, data_palette, render_type,
                default_color)

        self.__render(seq, ss, final_color_rbgs, filename, self.__draw_params)

        cairosvg.svg2png(url=filename+".svg", write_to=filename+".png",
                         output_width=self.__draw_params.output_width,
                         output_height=self.__draw_params.output_height)
        os.remove(filename+".svg")


    def __setup(self, draw_params):
        if draw_params is not None:
            self.__draw_params = draw_params


    def __get_rgb_colors(self, ss, seq, color_str, data, data_palette,
                         render_type, default_color):
        if default_color is None:
            default_color = COLORS["e"]
        return self.__colorer.get_rgb_colors(
                seq, ss, color_str, data, data_palette, render_type, default_color)


    def __render(self, seq, ss, colors, filename, params):
        r = render_rna.RNARenderer()

        pairmap = render_rna.get_pairmap_from_secstruct(ss)
        pairs = []
        for i in range(len(pairmap)):
            if pairmap[i] > i:
                pairs.append({"from": i, "to": pairmap[i], "p": 1.0,
                              "color": COLORS["e"]})

        r.setup_tree(ss, params.NODE_R, params.PRIMARY_SPACE, params.PAIR_SPACE)
        size = r.get_size()

        cell_size = max(size) + params.CELL_PADDING * 2

        svgobj = svg.svg("%s.svg" % filename, cell_size, cell_size)
        r.draw(svgobj, params.CELL_PADDING, params.CELL_PADDING, colors,
               pairs, seq, params.RENDER_IN_LETTERS)



def rna_draw(ss, seq=None, filename='secstruct', color_str=None, render_type=None,
             data=None, data_palette=None, data_file=None, default_color=None):

    rd = RNADrawer()
    render_type = get_render_type(render_type)

    # TODO finish data file
    if data_file is not None:
        df = pd.read_csv(data_file)
        data = np.zeros(len(seq))

    if data_palette is not None:
        data_palette = matplotlib.cm.get_cmap(data_palette)

    rd.draw(ss, seq, filename, color_str, render_type, data, data_palette,
            default_color)



def get_render_type(render_type_name):
    if render_type_name is None:
        return None

    render_type_name = render_type_name.lower()
    if   render_type_name == "res_type":
        return colorer.RenderType.RES_TYPE
    elif render_type_name == "paired":
        return colorer.RenderType.PAIRED
    elif render_type_name == "none":
        return None
    else:
        raise ValueError("unknown render type: " + render_type_name)


def main():
    args = parse_args()
    rd = RNADrawer()
    #rd.draw()



if __name__ == "__main__":
    main()



















