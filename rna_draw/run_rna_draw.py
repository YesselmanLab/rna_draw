import sys
import re
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
    ss, seq = parse_dbn_file(dbn_file_path)
    print('dbn path', dbn_file_path)
    print(ss)
    print(seq)
    overlap_count = rd.rna_draw(ss=ss, seq=seq, render_type='res_type', cluster=True)
    print("Successfully Drew Structure with Response Overlap:", overlap_count)

if __name__ == "__main__":
    dbn_file_path = sys.argv[1]
    run_rna_draw(dbn_file_path)
