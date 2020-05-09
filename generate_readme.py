from rna_draw import render, colorer

import markdown_generator as mg
import pandas as pd




def no_render_type_draw(rr, row, writer, i):
    rr.render(row["structure"], row["sequence"], row["colors"],
                          filename="resources/imgs/test_{i}".format(i=i))
    img = mg.Image(("resources/imgs/test_{i}.png".format(i=i)), "")
    writer.writeline(img)


def main():
    f = open("resources/md_files/install.md")
    lines = f.readlines()
    f.close()

    df = pd.read_csv("inputs.csv")
    rr = render.RNARender()

    with open('README.md', 'w') as f:
        writer = mg.Writer(f)
        writer.writelines(lines)
        for i, row in df.iterrows():
            if pd.isna(row["render_type"]):
                no_render_type_draw(rr, row, writer, i)
                continue

            if pd.notna(row["data"]):
                data = [float(x) for x in row["data"].split(";")]
                data_colors = colorer.color_by_data(data)
                print(data_colors)
                rr.render(row["structure"], row["sequence"], colors=data_colors,
                          filename="resources/imgs/test_{i}".format(i=i))

            else:
                rr.render(row["structure"], row["sequence"],
                            filename="resources/imgs/test_{i}".format(i=i))

            img = mg.Image(("resources/imgs/test_{i}.png".format(i=i)), "")
            writer.writeline(img)


if __name__ == "__main__":
    main()
