import rna_draw as rd

import markdown_generator as mg
import pandas as pd
import numpy as np


def get_draw_rna_cmd(args):
    s = "rna_draw "
    for a in args:
        s += "-{} {} ".format(a[0], a[1])
    return s


def get_draw_rna_py_cmd(args):
    s = "rd.rna_draw("
    for i, a in enumerate(args):
        s += "{}={}".format(a[0], a[1])
        if i != len(args)-1:
            s += ", "
    s += ")"
    return s

def get_args_from_dict(d):
    data = []
    for k, v in d.items():
        if v is None:
            continue
        if k == "ss":
            v = '"'+v+'"'
        data.append([k, v])
    return data

def get_args_from_dict_py(d):
    data = []
    for k, v in d.items():
        if v is None:
            continue
        if k == "data_vmin" or k == "data_vmax":
            pass
        else:
            v = "'"+v+"'"

        data.append([k, v])
    return data


def get_args_from_file(fname):
    f = open(fname)
    lines = f.readlines()
    f.close()

    args = {
        'ss'  : None,
        'seq' : None,
        'color_str' : None,
        'render_type' : None,
        'default_color' : None,
        'data_str' : None,
        'data_file' : None,
        'data_palette' : None,
        'data_vmin' : None,
        'data_vmax' : None,
        'data_ignore_restype' : None
    }
    for l in lines:
        spl = l.split("=")
        args[spl[0]] = spl[1].rstrip()
    return args


def underscore_to_spaces(s):
    spl = s.split("_")
    return " ".join(spl)


def render_header_title(name, level):
    tabs = "\t"*level
    spl = name.split("_")
    text_name = " ".join(spl)
    link_name = "#" + "-".join(spl)
    return tabs + "- [{}]({})\n".format(text_name, link_name)


def main():
    f = open("resources/md_files/install.md")
    lines = f.readlines()
    f.close()


    f = open("README.md", 'w')
    f.writelines(lines)
    writer = mg.Writer(f)

    df = pd.read_csv("readme_examples.csv")
    header_lines = []
    example_type = ""
    for i, row in df.iterrows():
        if row["example_type"] != example_type:
            header_lines.append(render_header_title(row["example_type"], 0))
            example_type = row["example_type"]

        header_lines.append(render_header_title(row["name"], 1))

    f.writelines(header_lines)

    example_type = ""
    for i, row in df.iterrows():
        if row["example_type"] != example_type:
            writer.write_heading(underscore_to_spaces(row["example_type"]), 2)
            example_type = row["example_type"]

        writer.write_heading(underscore_to_spaces(row["name"]), 3)
        dict_vals = get_args_from_file(
                "examples/" + row["example_type"] + "/" + row["name"] + ".txt")
        comment = dict_vals["comment"]
        del dict_vals["comment"]

        writer.writeline(comment)
        code = mg.Code("shell")
        cmd_args = get_args_from_dict(dict_vals)
        code.append(get_draw_rna_cmd(cmd_args))
        writer.write(code)
        code_py = mg.Code("python")
        code_py.append(get_draw_rna_py_cmd(get_args_from_dict_py(dict_vals)))
        writer.write(code_py)
        rd.rna_draw(
                    out="resources/imgs/test_{}".format(i),
                    **dict_vals)
        img = mg.Image(("resources/imgs/test_{}.png".format(i)), "")
        writer.writeline(img)
    #writer.writeline("[TOP](#how-to)")





if __name__ == "__main__":
    main()
