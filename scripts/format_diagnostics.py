import argparse
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-file", required=True, type=str, help="Input file")
args = parser.parse_args()

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

# all seqs should be of same length
assert(len(set([len(x[HEADER_TO_COL_IDX["sequence"]]) for x in d[1:]]))==1)
SEQLEN = len(d[1][HEADER_TO_COL_IDX["sequence"]])

# coords should be distinct if provided
all_coords = [x[HEADER_TO_COL_IDX['coord']] for x in d[1:]]
if len(all_coords) != len(set(all_coords)):
    print("Distinct coords : {}".format(len(set(all_coords))))

coord_regex= re.compile("chr[0-9XYM]+:[0-9]+-[0-9]+")
pos_regex= re.compile("[0-9]+-[0-9]+")

gc_counts = []
cg_counts = []

for j in range(1, len(d)):
    if any([d[j][HEADER_TO_COL_IDX[x]]=='' for x in ['motif_positions', 'motif_names', 'motif_orientations']]):
        if not all([d[j][HEADER_TO_COL_IDX[x]]=='' for x in ['motif_positions', 'motif_names', 'motif_orientations']]):
            print(j, d[j])
            print("If either of positions/names/orientations is empty, all must be empty, i.e. no motif introduced")
            exit(1)
        else:
            # can skip these lines, presumably background
            continue

    num_pos = len(d[j][HEADER_TO_COL_IDX["motif_positions"]].split(';'))
    num_names = len(d[j][HEADER_TO_COL_IDX["motif_names"]].split(';'))
    num_or = len(d[j][HEADER_TO_COL_IDX["motif_orientations"]].split(';'))

    if not all([pos_regex.match(x) for x in d[j][HEADER_TO_COL_IDX["motif_positions"]].split(';')]):
        print(j, d[j])
        print("Positions incorrectly formatted, should be start1-end1;start2-end2...if multiple else start-end")
        exit(1)

    if not (num_pos == num_names == num_or):
        print(j, d[j])
        print("Number of positions, number of names and number of orientations not matching!")
        exit(1)

    if not coord_regex.match(d[j][HEADER_TO_COL_IDX["coord"]]):
        print(j, d[j])
        print("Coord misspecified. Should be of form chr:start-end")
        exit(1)

    # check if there's a GC close enough to motif (if there is at least one motif)
    pos = [[int(z) for z in y.split("-")] for y in d[j][HEADER_TO_COL_IDX["motif_positions"]].split(";")]
    ors = d[j][HEADER_TO_COL_IDX["motif_orientations"]].split(";")

    for i, y in enumerate(pos):
        seq = d[j][HEADER_TO_COL_IDX["sequence"]][y[0]-5:y[1]+5]
        if "GC" not in seq:
            print(j, d[j], y[0], y[1], seq)
            print("GC missing near motif instance")
            exit(1)
        

    gc_counts.append(d[j][HEADER_TO_COL_IDX["sequence"]].upper().count("GC"))
    cg_counts.append(d[j][HEADER_TO_COL_IDX["sequence"]].upper().count("CG"))

gc_counts = np.array(gc_counts)
cg_counts = np.array(cg_counts)

print("GC% Mean : {:.2f} Std : {:.2f}".format(100*np.mean(gc_counts/SEQLEN), np.std(100*gc_counts/SEQLEN)))
print("CG% Mean : {:.2f} Std : {:.2f}".format(100*np.mean(cg_counts/SEQLEN), np.std(100*cg_counts/SEQLEN)))
