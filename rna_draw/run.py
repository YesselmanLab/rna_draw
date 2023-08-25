#import rna_draw as rd

###rd.rna_draw(
    #ss="(((.(((..((((........))))..((((.......))))......(((((.......)))))))).)))....",
    #seq="GCGGAUUUAGCUCAGCCGGGAGAGCGCCAGACUGGACACCUGGAGGUCCUGUGCCCGAUCCACAGAAUUCGCACCA",
    #render_type="res_type",
#)
# rd.rna_draw(ss=".(((.....))).")

import os
import re
from PIL import Image
import rna_draw as rd

# Location of the DBN files (relative path from the Desktop)
DBN_FOLDER = "dbnFiles"

# Function to create an empty image
def create_empty_image(filename):
    img = Image.new('RGB', (800, 600), color='white')
    img.save(filename)

# Iterate through each .dbn file in the folder
for dbn_file_name in os.listdir(DBN_FOLDER):
    if dbn_file_name.endswith(".dbn"):
        dbn_file_path = os.path.join(DBN_FOLDER, dbn_file_name)

        # Extract the name of the DBN file without the extension
        file_name_no_ext = os.path.splitext(dbn_file_name)[0]

        try:
            with open(dbn_file_path, "r") as file:
                content = file.read()

                # Extract the sequence of nucleotides
                seq = re.search(r"^[ACGUI]+$", content, re.MULTILINE)
                seq = seq.group(0) if seq else ""

                # Replace 'I' with 'A' in the sequence
                seq = seq.replace('I', 'A') # Disable this line if needed

                # Extract the sequence of parenthesis, dots, and square brackets
                ss = re.search(r"^[().\[\]]+$", content, re.MULTILINE)
                ss = ss.group(0) if ss else ""

                # Replace square brackets with dots
                ss = ss.replace('[', '.').replace(']', '.')

                # Check if the sequences contain only valid characters
                if not seq or not ss:
                    raise ValueError("Invalid characters detected.")

                # Run the rna_draw command with the extracted sequences
                rd.rna_draw(ss=ss, seq=seq, render_type="res_type")

            # Copy the image to the destination folder with the same name as the DBN file
            os.rename("secstruct.png", f"Images/{file_name_no_ext}.png")
        except Exception as e:
            print(f"Error processing {dbn_file_name}: {str(e)}. Creating empty image.")
            create_empty_image(f"Images/{file_name_no_ext}-ERROR.png")



print("Processing complete.")

