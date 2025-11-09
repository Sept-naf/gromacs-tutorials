#!/usr/bin/env python3
"""
water_deletor.py  –  remove waters inside membrane slab & write new topology
default files:  -in solv.gro  -out solv_delW.gro  -top topol.top  -topout topol_delW.top
usage:  python water_deletor.py  [-ref N] [-middle N] [-nwater 3] [-v]
"""
import sys, shutil, argparse
import re


parser = argparse.ArgumentParser(description="Delete waters inside membrane slab & write new topology")
parser.add_argument("-in",      dest="ingro",   default="solv.gro",          help="input .gro file  (default: solv.gro)")
parser.add_argument("-out",     dest="outgro",  default="solv_delW.gro",     help="output .gro file (default: solv_delW.gro)")
parser.add_argument("-top",     dest="intop",   default="topol.top",         help="input topol.top  (default: topol.top)")
parser.add_argument("-topout",  dest="outtop",  default="topol_delW.top",    help="output topology  (default: topol_delW.top)")
parser.add_argument("-ref",     dest="refatom", default="N",                 help="atom name for bilayer top/bottom (default: N)")
parser.add_argument("-middle",  dest="midatom", default="N",                 help="atom name for bilayer mid-plane (default: N)")
parser.add_argument("-nwater",  dest="nwater",  type=int, default=3,         help="atoms per water (default: 3)")
parser.add_argument("-v",       dest="verb",    action="store_true",         help="verbose: print deleted residues")
args = parser.parse_args()

# -------------------------------------------------
# 1  read .gro
# -------------------------------------------------
with open(args.ingro) as fh:
    lines = fh.readlines()

header      = lines[0]
natoms_line = lines[1]
boxline     = lines[-1]
atomlines   = lines[2:-1]

# -------------------------------------------------
# 2  helpers  –  WORK IN NM (gro file unit)
# -------------------------------------------------
def get_z(l):
    """z coordinate in nm from .gro line (cols 36-44)"""
    return float(l[36:44])

# -------------------------------------------------
# 3  find POPC phosphate N on both leaflets
# -------------------------------------------------
popc_lines = [l for l in atomlines if l[5:10].strip() == "POPC"]
n_lines    = [l for l in popc_lines if l[12:16].strip() == "N"]

if not n_lines:
    sys.exit("No POPC phosphate nitrogen (atom name N) found – check residue/atom names")

z_mid = sum(get_z(l) for l in n_lines) / len(n_lines)
print(f" Bilayer middle (avg POPC N) z = {z_mid:.3f} nm")

top_n = [l for l in n_lines if get_z(l) > z_mid]
bot_n = [l for l in n_lines if get_z(l) < z_mid]

if not top_n or not bot_n:
    sys.exit("Could not split POPC N into top/bottom – check bilayer sanity")

z_top_avg = sum(get_z(l) for l in top_n) / len(top_n)
z_bot_avg = sum(get_z(l) for l in bot_n) / len(bot_n)

print(f" Top POPC N avg z = {z_top_avg:.3f} nm")
print(f" Bottom POPC N avg z = {z_bot_avg:.3f} nm")

# -------------------------------------------------
# 4  scan waters & mark deletions
# -------------------------------------------------

ow_indices = [i for i, l in enumerate(atomlines)
              if l[5:10].strip() == "SOL" and re.search(r"\bOW\b", l[10:15])]
#print(ow_indices)


to_delete = set()
for idx in ow_indices:
    z = get_z(atomlines[idx])
    #print(z)
    if z_bot_avg < z < z_top_avg:
        for j in range(args.nwater):
            to_delete.add(idx + j)

n_wat_del = len(to_delete) // args.nwater
print(f" {n_wat_del} water molecules inside membrane (between POPC N layers) – deleting")

# -------------------------------------------------
# 5  build new .gro list & count surviving SOL molecules
# -------------------------------------------------
new_lines = [l for i, l in enumerate(atomlines) if i not in to_delete]

# --- count SOL residues ---
sol_set = set()
for l in new_lines:
    if l[5:10].strip() == "SOL":
        sol_set.add(l[0:15])          # unique residue key: resnum+resname+atomname
new_sol = len(sol_set) // args.nwater   # whole waters

new_natoms = len(new_lines)

# -------------------------------------------------
# 6  write new .gro
# -------------------------------------------------
with open(args.outgro, "w") as fh:
    fh.write(header)
    fh.write(f"{new_natoms:>10d}\n")
    for l in new_lines:
        fh.write(l)
    fh.write(boxline)
print(f" Output .gro written to {args.outgro}")

# -------------------------------------------------
# 7  patch topology – keep only the *first* SOL line
# -------------------------------------------------
shutil.copy(args.intop, args.outtop)

inside_molecules = False
sol_patched = False
with open(args.outtop) as fh, open(args.outtop + ".tmp", "w") as out:
    for line in fh:
        if line.strip().startswith("[ molecules ]"):
            inside_molecules = True
            out.write(line)
            continue
        # skip *any* subsequent SOL line
        if inside_molecules and line.strip().startswith("SOL"):
            if not sol_patched:               # first SOL → update number
                parts = line.split()
                out.write(f"{parts[0]:<10} {new_sol}\n")
                sol_patched = True
            # **all other SOL lines are ignored (deleted)**
            continue
        # once we leave the section, stop filtering
        if inside_molecules and line.strip().startswith("[") and not line.strip().startswith("[ molecules ]"):
            inside_molecules = False
        out.write(line)

shutil.move(args.outtop + ".tmp", args.outtop)
print(f" Patched {args.outtop}  –  SOL count = {new_sol}  (extra SOL lines removed)")
print(" Done.")
