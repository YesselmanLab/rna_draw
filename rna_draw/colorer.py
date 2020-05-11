import matplotlib.colors
import matplotlib.cm

import numpy as np
import seaborn as sns

COLORS = {"r": [255, 102, 102],
          "g": [113, 188, 120],
          "b": [51, 153, 255],
          "k": [1, 0, 0],
          "y": [255, 211, 0],
          "c": [0, 255, 255],
          "m": [255, 0, 255],
          "w": [255, 255, 255],
          "e": [150, 150, 150],
          "o": [231, 115, 0],
          "i": [51, 204, 204],
          "h": [51, 153, 255]}

class RenderType:
    RES_TYPE = 0
    PAIRED = 1
    MOTIF = 2


class Colorer(object):
    def __init__(self):
        self.seq = None
        self.ss = None
        self.color_str = None
        self.data = None
        self.render_type = None


    def get_rgb_colors(self, seq, ss, color_str=None, data=None, data_palette=None,
                       render_type=None, default_color=COLORS["e"]):
        """
        :param seq: sequence of the RNA
        :param ss:  dot bracket notation for sec
        :param color_str: string to store color names, or color ranges
        :param data: data
        :param data_palette:
        :param render_type:
        :param default_color:
        :return:
        """

        self.seq = seq
        self.ss = ss
        self.color_str = color_str
        self.data = data
        self.render_type = render_type
        self.set_colors = np.zeros(len(seq))
        self.rgb_colors = [default_color for i in range(len(seq))]

        if len(seq) != len(ss):
            raise ValueError("sequence and structure must be the same length")

        if render_type is not None:
            rgb_colors = []
            set_colors = np.ones(len(seq))
            if   render_type == RenderType.RES_TYPE:
                rgb_colors = self.__color_by_restype()
            elif render_type == RenderType.PAIRED:
                rgb_colors = self.__color_by_paired()
            self.__add_rgb_colors(rgb_colors, set_colors)

        if data is not None:
            rgb_colors = color_by_data(data, data_palette)
            set_colors = np.ones(len(seq))
            self.__add_rgb_colors(rgb_colors, set_colors)

        if color_str is not None:
            rgb_colors, set_colors = self.__parse_color_str(color_str)
            self.__add_rgb_colors(rgb_colors, set_colors, override=True)

        return self.rgb_colors


    def __parse_color_str(self, color_str):
        rgb_colors = [COLORS["e"] for i in range(len(self.seq))]

        spl = color_str.split(";")
        contains_digit = any(map(str.isdigit, color_str))
        # there is a single color code for each sequence position
        if len(spl) == 1 and not contains_digit:
            if len(color_str) != len(self.seq):
                raise ValueError(
                        "there are no ; in color str thus it must match the seq length")
            return parse_color_single_letter_codes(color_str), np.ones(len(self.seq))

        elif not contains_digit:
            if len(spl) != len(self.seq):
                raise ValueError(
                        "no residue numbers are specified in color str, must match sequence length")
        set_colors = np.zeros(len(self.seq))

        if len(spl) == len(self.seq):
            for i, color_code in enumerate(spl):
                rgb_colors[i] = self.__parse_color_code(color_code)
                set_colors[i] = 1
            return rgb_colors, set_colors

        for color_code_range in spl:
            if len(color_code_range) < 2:
                continue
            num_range, color_code = color_code_range.split(":")
            rgb_color = self.__parse_color_code(color_code)
            if len(num_range.split("-")) > 1:
                min_num, max_num = [int(x) for x in num_range.split("-")]
                if min_num > max_num:
                    raise ValueError(
                            "when supplying a range of res nums for colors smaller number "+
                            "must come first")
            else:
                min_num, max_num = int(num_range), int(num_range)

            for i in range(min_num-1, max_num):
                if set_colors[i]:
                    raise ValueError("position {} has two colors assigned to it".format(i))
                rgb_colors[i] = rgb_color
                set_colors[i] = 1

        return rgb_colors, set_colors


    def __parse_color_code(self, color_code):
        if len(color_code) == 1:
            return COLORS[color_code]

        elif color_code in sns.xkcd_rgb:
            return xkcd_color_name_to_rgb(color_code)

        else:
            raise ValueError("color_code: " + color_code + " is unknown")


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
        rgb_colors = []
        for e in self.seq:
            if   e == "A":
                rgb_colors.append(COLORS["y"])
            elif e == "C":
                rgb_colors.append(COLORS["g"])
            elif e == "G":
                rgb_colors.append(COLORS["r"])
            elif e == "U" or e == "T":
                rgb_colors.append(COLORS["b"])
            elif e == "N":
                rgb_colors.append(COLORS["e"])
            else:
                rgb_colors.append(COLORS["m"])
        return rgb_colors


    def __color_by_paired(self):
        rgb_colors = []
        for e in self.ss:
            if e == '.':
                rgb_colors.append(COLORS["y"])
            else:
                rgb_colors.append(COLORS["b"])
        return rgb_colors


def parse_color_single_letter_codes(color_string):
    colors = []
    for e in color_string:
        colors.append(COLORS[e])
    return colors


def color_by_res_type(seq, ss):
    colors = []
    for e in seq:
        if e == "A":
            colors.append(COLORS["y"])
        elif e == "C":
            colors.append(COLORS["g"])
        elif e == "G":
            colors.append(COLORS["r"])
        elif e == "U" or e == "T":
            colors.append(COLORS["b"])
        else:
            colors.append(COLORS["e"])
    return colors


def color_by_data(data, palette=None):
    if palette is None:
        palette = matplotlib.cm.get_cmap("Reds")

    norm = matplotlib.colors.Normalize(vmin=min(data), vmax=max(data))
    cmap = palette

    colors = []
    for d in data:
        rgb = [int(255*x) for x in cmap(norm(d))[0:3]]

        colors.append(rgb)
        #colors.append(matplotlib.colors.to_rgb(cmap(norm(d))))
    return colors


def xkcd_color_name_to_rgb(name):
    raw_rgb = matplotlib.colors.hex2color(sns.xkcd_rgb[name])
    return [int(255*x) for x in raw_rgb]

