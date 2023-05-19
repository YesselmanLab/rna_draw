import matplotlib.colors as mcolors
import re
import rle
from pprint import pprint
import numpy as np
import seaborn as sns

from rna_draw.parameters import RenderType
from rna_draw.data import Data

DEFAULT_COLORS = {
    "r": [255 / 255, 102 / 255, 102 / 255],
    "g": [113 / 255, 188 / 255, 120 / 255],
    "b": [51 / 255, 153 / 255, 255 / 255],
    "y": [255 / 255, 211 / 255, 0 / 255],
    "c": [0 / 255, 255 / 255, 255 / 255],
    "m": [255 / 255, 0 / 255, 255 / 255],
    "w": [255 / 255, 255 / 255, 255 / 255],
    "e": [150 / 255, 150 / 255, 150 / 255],
    "o": [231 / 255, 115 / 255, 0 / 255],
}



class Colorer(object):

    COLORS = DEFAULT_COLORS

    def __init__(self):
        self.seq = None
        self.ss = None
        self.color_str = None
        self.data = None
        self.render_type = None
        self.palette = self.__class__.COLORS

    def get_rgb_colors(
        self, seq, ss, color_str=None, data=None, render_type=None, default_color=None
    ):
        if default_color is None:
            default_color = DEFAULT_COLORS["e"]

        if len(seq) != len(ss):
            raise ValueError("sequence and structure must be the same length")

        self.seq = seq
        self.ss = ss
        self.color_str = color_str
        self.data = data
        self.render_type = render_type
        self.set_colors = np.zeros(len(seq))
        self.rgb_colors = [default_color for i in range(len(seq))]
        self.default_color = default_color

        if render_type is not None and data is not None:
            raise ValueError("cannot set both render_type and data for the same res")

        if render_type is not None:
            rgb_colors = []
            set_colors = np.ones(len(ss))
            if render_type == RenderType.RES_TYPE:
                rgb_colors = self.__color_by_restype()
            elif render_type == RenderType.PAIRED:
                rgb_colors = self.__color_by_paired()
            elif render_type == RenderType.STRAND:
                rgb_colors = self.__color_by_strand()
            self.__add_rgb_colors(rgb_colors, set_colors)

        if data is not None:
            rgb_colors = color_by_data(data)
            set_colors = np.ones(len(ss))

            for i, e in enumerate(self.seq):
                if e in data.ignore_restype:
                    rgb_colors[i] = self.default_color

            self.__add_rgb_colors(rgb_colors, set_colors, override=True)

        if color_str is not None:
            rgb_colors, set_colors = self.__parse_color_str(color_str)
            self.__add_rgb_colors(rgb_colors, set_colors, override=True)

        return self.rgb_colors

    @classmethod
    def color_palette(cls,palette='colorblind'):
        VALID_PALETTES = {"deep", "muted", "pastel", "bright", "dark", "colorblind"}
        if palette not in VALID_PALETTES:
            raise ValueError("color_palette: palette must be one of %r." % VALID_PALETTES)
        cls.COLORS = dict(zip(['b', 'o', 'g', 'r', 'p', 'R', 'm', 'e', 'y', 'c'],
                                 list(map(mcolors.hex2color,sns.color_palette(palette,as_cmap=True)))))

    def color_palette(self,palette='colorblind'):
        VALID_PALETTES = {"deep", "muted", "pastel", "bright", "dark", "colorblind"}
        if palette not in VALID_PALETTES:
            raise ValueError("color_palette: palette must be one of %r." % palette)
        self.palette = dict(zip(['b', 'o', 'g', 'r', 'p', 'R', 'm', 'e', 'y', 'c'],
                                list(map(mcolors.hex2color,sns.color_palette(palette,as_cmap=True)))))
        
    def __parse_color_str(self, color_str):
        rgb_colors = [DEFAULT_COLORS["e"] for i in range(len(self.seq))]

        spl = color_str.split(";")
        contains_digit = any(map(str.isdigit, color_str))
        # there is a single color code for each sequence position
        if len(spl) == 1 and not contains_digit:
            if len(color_str) != len(self.seq):
                raise ValueError(
                    "there are no ; in color str thus it must match the secondary structure length"
                )
            return parse_color_single_letter_codes(color_str), np.ones(len(self.ss))

        elif not contains_digit:
            if len(spl) != len(self.ss):
                raise ValueError(
                    "no residue numbers are specified in color str, must match secondary structure length"
                )
        set_colors = np.zeros(len(self.ss))

        # there is only one color per residue
        if len(spl) == len(self.ss):
            for i, color_code in enumerate(spl):
                rgb_colors[i] = parse_color_code(color_code)
                set_colors[i] = 1
            return rgb_colors, set_colors

        for color_code_range in spl:
            if len(color_code_range) < 2:
                continue
            num_range, color_code = color_code_range.split(":")
            rgb_color = parse_color_code(color_code)
            if len(num_range.split("-")) > 1:
                min_num, max_num = [int(x) for x in num_range.split("-")]
                if min_num > max_num:
                    raise ValueError(
                        "when supplying a range of res nums for colors smaller number "
                        + "must come first"
                    )
            else:
                min_num, max_num = int(num_range), int(num_range)

            for i in range(min_num - 1, max_num):
                if set_colors[i]:
                    raise ValueError(
                        "position {} has two colors assigned to it".format(i)
                    )
                rgb_colors[i] = rgb_color
                set_colors[i] = 1

        return rgb_colors, set_colors

    def __add_rgb_colors(self, rgb_colors, set_colors, override=False):
        for i, rgb_color in enumerate(rgb_colors):
            if not set_colors[i]:
                continue
            if self.set_colors[i] and not override:
                raise ValueError("cannot set color it would override existing color")
            else:
                self.rgb_colors[i] = rgb_color
                self.set_colors[i] = 1

    def __color_by_restype(self):
        res_colors = {
            "A" : "y",
            "C" : "g",
            "G" : "r",
            "U" : "b",
            "T" : "b",
            "N" : "e"}
        return map(lambda s: self.palette[res_colors.get(s,"m")],self.seq)


    def __color_by_paired(self):
        return map(lambda ss: self.palette["y" if ss == "." else "b"],self.ss)

    def __color_by_strand(self):
        colors = list(self.palette.keys())
        rgb_colors = [self.palette[colors[0]] for i in range(len(self.ss))]
        breaks = {x.start() for x in re.finditer('\+', self.ss)}
        for index,element in enumerate(breaks):
            for i in range(element,len(rgb_colors)):
                rgb_colors[i] = self.palette[colors[index + 1]]
        return rgb_colors


def parse_color_single_letter_codes(color_string):
    colors = []
    for e in color_string:
        colors.append(self.palette[e])
    return colors


def parse_color_code(color_code):
    if len(color_code) == 1:
        return self.palette[color_code]

    elif color_code in sns.xkcd_rgb:
        return xkcd_color_name_to_rgb(color_code)

    else:
        raise ValueError("color_code: " + color_code + " is unknown")


def color_by_data(data):
    norm = matplotlib.colors.Normalize(vmin=data.min, vmax=data.max)
    colors = []
    for i, d in enumerate(data.vals):
        rgb = [x for x in data.palette(norm(d))[0:3]]
        colors.append(rgb)
    return colors


def xkcd_color_name_to_rgb(name):
    raw_rgb = matplotlib.colors.hex2color(sns.xkcd_rgb[name])
    return raw_rgb


