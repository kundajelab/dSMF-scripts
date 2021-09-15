import argparse

# extracts and prints all motif instances 

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-file", required=True, type=str, help="Input file")
parser.add_argument("-r", "--rev-comp", default=1, type=int, help="Return same strand of motif instances")
args = parser.parse_args()

def rc(s):
    RC={'A':'T','C':'G','G':'C','T':'A'}
    return ''.join(RC[x] for x in s)[::-1]

REQD_HEADERS = ["sequence", "motif_positions", "motif_names", "motif_orientations", "coord"]

HEADER_TO_COL_IDX = dict()

with open(args.input_file) as f:
    d = [x.strip().split(',') for x in f]

d = [[y.strip() for y in x] for x in d]

for h in REQD_HEADERS:
    if h not in d[0]:
        print("Header {} required!".format(h))
        exit(1)
    HEADER_TO_COL_IDX[h] = d[0].index(h)

for x in d[1:]:
    pos = [[int(z) for z in y.split("-")] for y in x[HEADER_TO_COL_IDX["motif_positions"]].split(";") if y!='']
    ors = x[HEADER_TO_COL_IDX["motif_orientations"]].split(";")

    for i, y in enumerate(pos):
        seq = x[HEADER_TO_COL_IDX["sequence"]][y[0]:y[1]]

        if args.rev_comp and ors[i]=="-":
            seq = rc(seq)
        print(seq)

