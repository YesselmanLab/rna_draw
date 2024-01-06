import os
import re
import traceback
from PIL import Image
from rna_draw import draw
rna_draw = draw.RNADrawer()

DBN_FOLDER = "dbnFiles"
process_all_files = False
specific_file_name = "bpRNA_CRW_15130.dbn"

def create_empty_image(filename):
    img = Image.new('RGB', (800, 600), color='white')
    img.save(filename)

def process_file(dbn_file_name):
    dbn_file_path = os.path.join(DBN_FOLDER, dbn_file_name)
    file_name_no_ext = os.path.splitext(dbn_file_name)[0]

    try:
        with open(dbn_file_path, "r") as file:
            content = file.read()
            seq = re.search(r"^[ACGUI]+$", content, re.MULTILINE)
            seq = seq.group(0) if seq else ""
            seq = seq.replace('I', 'A')
            ss = re.search(r"^[().\[\]]+$", content, re.MULTILINE)
            ss = ss.group(0) if ss else ""
            ss = ss.replace('[', '.').replace(']', '.').replace('}', '.').replace('{', '.').replace('<', '.').replace('>', '.')
            
            if not seq or not ss:
                raise ValueError("Invalid characters detected.")

            rna_draw.draw(ss=ss, seq=seq, render_type='res_type')

        os.rename("secstruct.png", f"Images/{file_name_no_ext}.png")
    except SystemExit:
        print(f"A sys.exit() call was caught while processing {dbn_file_name}. Creating empty image.")
        create_empty_image(f"Images/{file_name_no_ext}-ERROR.png")
    except Exception as e:
        tb = traceback.extract_tb(e.__traceback__)
        filename, line, func, text = tb[-1]
        print(f"Error processing {dbn_file_name}: {str(e)}. Creating empty image. Error at line {line}")
        create_empty_image(f"Images/{file_name_no_ext}-ERROR.png")

if process_all_files:
    for dbn_file_name in os.listdir(DBN_FOLDER):
        if dbn_file_name.endswith(".dbn"):
            process_file(dbn_file_name)
else:
    if specific_file_name.endswith(".dbn"):
        process_file(specific_file_name)

print("Processing complete.")
