import sys
import re
import traceback
import rna_draw as rd

def parse_dbn_file(dbn_file_path):
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

    return ss, seq

def run_rna_draw(dbn_file_path):
    try:
        ss, seq = parse_dbn_file(dbn_file_path)
        overlap_count = rd.rna_draw(ss=ss, seq=seq, render_type='res_type', cluster='True', out=dbn_file_path.replace('.dbn', ''))
    except Exception as e:
        error_message = f"An error occurred: {e}\n"
        error_traceback = traceback.format_exc()
        
        error_file_path = '/work/yesselmanlab/nklein/Fails/' + dbn_file_path.replace('.dbn', '_error.txt')

        with open(error_file_path, 'w') as error_file:
            error_file.write(error_message)
            error_file.write(error_traceback)
    except SystemExit:
        error_message = f"An error occurred: System Exit\n"
        
        error_file_path = '/work/yesselmanlab/nklein/Fails/' + dbn_file_path.replace('.dbn', '_error.txt')

        with open(error_file_path, 'w') as error_file:
            error_file.write(error_message)

if __name__ == "__main__":
    dbn_file_path = sys.argv[1]
    run_rna_draw(dbn_file_path)
