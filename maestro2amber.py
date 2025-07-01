import sys
from collections import defaultdict

inputpdb = sys.argv[1]
outputpdb = sys.argv[2]
fr = open(inputpdb, "r")
fw = open(outputpdb, "w")

# Read all lines and group ATOM lines by residue
residues = defaultdict(list)
lines = []
for line in fr:
    lines.append(line)
    if line.startswith('ATOM'):
        chain = line[21:22].strip()
        resnum = line[22:26].strip()
        inscode = line[26:27].strip()
        key = (chain, resnum, inscode)
        residues[key].append(line)

# Process each residue group and write lines back with modified residue names
i = 0
while i < len(lines):
    line = lines[i]
    if line.startswith('ATOM'):
        atom_name = line[12:16].strip()
        if atom_name.startswith('H'):
            i += 1
            continue  # Skip hydrogen atoms

        first_line_of_group = residues[(line[21:22].strip(), line[22:26].strip(), line[26:27].strip())][0]
        original_resname = first_line_of_group[17:20].strip()
        group_atom_names = [atom_line[12:16].strip() for atom_line in residues[(line[21:22].strip(), line[22:26].strip(), line[26:27].strip())]]
        new_resname = original_resname

        if original_resname == 'GLU' and 'HE2' in group_atom_names:
            new_resname = 'GLH'
        elif original_resname == 'ASP' and 'HD2' in group_atom_names:
            new_resname = 'ASH'
        elif original_resname == 'LYS' and 'HZ1' not in group_atom_names:
            new_resname = 'LYN'
        elif original_resname == 'CYS' and 'HG' not in group_atom_names:
            new_resname = 'CYX'
        elif original_resname == 'HIS':
            has_hd1 = 'HD1' in group_atom_names
            has_he2 = 'HE2' in group_atom_names
            if has_hd1 and has_he2:
                new_resname = 'HIP'
            elif has_hd1:
                new_resname = 'HID'
            elif has_he2:
                new_resname = 'HIE'

        modified_line = line[:17] + new_resname.ljust(3) + line[20:]
        fw.write(modified_line)

        # Process the rest of the atoms in the same residue
        while i+1 < len(lines) and lines[i+1].startswith('ATOM'):
            next_line = lines[i+1]
            next_atom_name = next_line[12:16].strip()
            if next_atom_name.startswith('H'):
                i += 1
                continue  # Skip hydrogen atoms

            if next_line[21:22].strip() == line[21:22].strip() and \
               next_line[22:26].strip() == line[22:26].strip() and \
               next_line[26:27].strip() == line[26:27].strip():
                modified_line = next_line[:17] + new_resname.ljust(3) + next_line[20:]
                fw.write(modified_line)
                i += 1
            else:
                break
    else:
        fw.write(line)
    i += 1

fr.close()
fw.close()
