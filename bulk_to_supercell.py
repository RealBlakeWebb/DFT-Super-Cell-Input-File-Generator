"""
Name: Blake Webb
Date: 6/26/2025
Project: Bulk to Super Cell input generator
Description: Program takes bulk input file and outputs scaled supercell input file
"""

from decimal import Decimal, getcontext
import os
import re

# Precision control for clean decimal formatting
getcontext().prec = 20

# ================================
# === HELPER FUNCTIONS ==========
# ================================
"""
ID: format_fractional_coord(val)
Parameters: 
- val - numerical value thats being reformated
Description: helper function to format reduced coordinates
"""
def format_fractional_coord(val):
    """Format a fractional or float string to 9 decimal places for QE."""
    if '/' in val:
        num, denom = map(int, val.split('/'))
        value = Decimal(num) / Decimal(denom)
    else:
        value = Decimal(str(val))
    return f"{value:.9f}"
"""
ID: read_input_file(filepath)
Parameters: 
- filepath
Description: helper function to open file.
"""
def read_input_file(filepath):
    with open(filepath, 'r') as f:
        return f.read()

"""
ID: extract_block(text, block_name)
Parameters: 
- text - file
- block_name - name of block
Description: helper function to autofill "CONTROL", "ELECTRONS", "IONS", "CELL".
"""
def extract_block(text, block_name):
    pattern = rf"&{block_name}(.*?)\/"
    match = re.search(pattern, text, re.DOTALL)
    return f"&{block_name}{match.group(1)}/\n" if match else ""
"""
FIX
"""
def extract_block(text, block_name):
    # Matches everything from &BLOCK_NAME to a line that only contains a slash
    pattern = rf"&{block_name}\b(.*?)(?:\n\s*/)"
    match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
    if match:
        return f"&{block_name}{match.group(1)}\n/"
    else:
        print(f"⚠️ Warning: Block &{block_name} not found.")
        return ""

"""
ID: extract_named_section(text, name)
Parameters: 
- text - file 
- name - name of block
Description: helper function to autofill "ATOMIC_POSITIONS (crystal)","K_POINTS {automatic}".
"""
def extract_named_section(text, name):
    pattern = rf"{name}.*?\n((?:.|\n)*?)(?:\n\n|\Z)"
    match = re.search(pattern, text, re.IGNORECASE)
    return f"{name}\n{match.group(1).strip()}" if match else ""
"""
ID: parse_atomic_positions(block)
Parameters: 
- block
Description: helper function to parse atomic positions from original input file.
"""
def parse_atomic_positions(block):
    lines = block.strip().splitlines()
    atoms = []
    for line in lines[1:]:  # skip ATOMIC_POSITIONS header
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms
"""
FIX
"""
def parse_atomic_positions(block):
    lines = block.strip().splitlines()
    atoms = []
    if not lines or not lines[0].lower().startswith("atomic_positions"):
        print("⚠️ Warning: ATOMIC_POSITIONS block not found or malformed.")
        return []
    for line in lines[1:]:
        if line.strip():
            parts = line.split()
            if len(parts) >= 4:
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            else:
                print(f"⚠️ Skipping malformed atomic line: {line}")
    return atoms

"""
ID: parse_celldm(block)
Parameters: 
- block
Description: helper function to parse atomic positions from original input file.
"""
def parse_celldm(block):
    celldm = {}
    matches = re.findall(r'celldm\((\d+)\)\s*=\s*([0-9.eE+-]+)',block)
    for idx, val in matches:
        celldm[int(idx)] = float(val)
    return celldm
"""
ID: scale_atomic_positions(atoms, sx, sy, sz)
Parameters: 
- atoms
- sx
- sy
- sz
Description: helper function to parse atomic positions from original input file.
"""
def scale_atomic_positions(atoms, sx, sy, sz):
    scaled_atoms = []
    for label, x, y, z in atoms:
        for i in range(sx):
            for j in range(sy):
                for k in range(sz):
                    new_x = (x + i) / sx
                    new_y = (y + j) / sy
                    new_z = (z + k) / sz
                    scaled_atoms.append((label, new_x, new_y, new_z))
    return scaled_atoms

"""
ID: format_atomic_positions(atoms)
Parameters: 
- atoms
Description: helper function to formatic new atomic positions for supercell input file.
"""
def format_atomic_positions(atoms):
    lines = ["ATOMIC_POSITIONS (crystal)"]
    for label, x, y, z in atoms:
        lines.append(f"{label:<3} {x:.9f} {y:.9f} {z:.9f}")
    return "\n".join(lines)
"""
ID: update_celldm(celldm1, celldm3, sx, sz)
Parameters: 
- celldm1, 
- celldm3, 
- sx
- sz
Description: helper function to update cell dimensions for supercell input file.
"""
def update_celldm(celldm1, celldm3, sx, sz):
    new_celldm1 = celldm1 * sx
    new_celldm3 = celldm3 * sz / sx
    return new_celldm1, new_celldm3
"""
ID: scale_kpoints(kpoints_line, sx, sy, sz)
Parameters: 
- kpoints_line
- sx
- sy
- sz
Description: helper function to update k points for supercell input file.
"""
def scale_kpoints(kpoints_line, sx, sy, sz):
    k = list(map(int, kpoints_line.split()))
    new_k = [max(1, k[0] // sx), max(1, k[1] // sy), max(1, k[2] // sz)] + k[3:]
    return f"K_POINTS {{automatic}}\n{' '.join(map(str, new_k))}"

# ================================
# === MAIN WORKFLOW =============
# ================================

def main():
    print("\n=== Quantum ESPRESSO Supercell Input Generator ===\n")

    # === Load and read original file ===
    input_file = input("Enter path to QE bulk input file: ").strip()
    text = read_input_file(input_file)

    # === Prompt for scaling and cleave direction ===
    scaling = input("Enter scaling [format: X,Y,Z]: ")
    sx, sy, sz = map(int, scaling.replace(",", " ").split())
    
    #when cleaving is later implemented
    #cleavedir = input("Enter cleave direction (x/y/z): ").lower()

    # === Extract relevant blocks from original ===
    control_block = extract_block(text, "CONTROL")
    electrons_block = extract_block(text, "ELECTRONS")
    ions_block = extract_block(text, "IONS")
    cell_block = extract_block(text, "CELL")
    atomic_species_block = extract_named_section(text, "ATOMIC_SPECIES")
    atomic_positions_block = extract_named_section(text, "ATOMIC_POSITIONS (crystal)")
    k_points_block = extract_named_section(text, "K_POINTS {automatic}")

    # === Parse and scale atomic positions ===
    atoms = parse_atomic_positions(atomic_positions_block)
    scaled_atoms = scale_atomic_positions(atoms, sx, sy, sz)
    nat = len(scaled_atoms)
    formatted_positions = format_atomic_positions(scaled_atoms)

    # === Extract and scale celldm(1) and celldm(3) from SYSTEM ===
    system_block = extract_block(text, "SYSTEM")
    celldm_dict = parse_celldm(system_block)
    celldm1_orig = celldm_dict.get(1, 1.0) #fills in 1.0 if value not found
    celldm3_orig = celldm_dict.get(3, 1.0) #fills in 1.0 if vlaue not found
    celldm1_new, celldm3_new = update_celldm(celldm1_orig, celldm3_orig, sx, sz)

    # === Get ntyp (number of atomic species) ===
    ntyp = sum(1 for line in atomic_species_block.splitlines() if line.strip() and not line.startswith("ATOMIC_SPECIES"))

    # === Format new SYSTEM block ===
    system_block_new = f"""&SYSTEM
  ibrav = 4,
  celldm(1) = {celldm1_new},
  celldm(3) = {celldm3_new},
  nat = {nat},
  ntyp = {ntyp},
  ecutwfc = 40,
  ecutrho = 320,
  occupations = 'smearing',
  smearing = 'gaussian',
  degauss = 0.015,
  input_dft = 'wc',
/
"""

    # === Adjust and format K_POINTS ===
    kpoints_line = k_points_block.strip().splitlines()[1]
    kpoints_block_new = scale_kpoints(kpoints_line, sx, sy, sz)

    # === Assemble all blocks ===
    output = "\n".join([
        control_block.strip(),
        system_block_new.strip(),
        electrons_block.strip(),
        ions_block.strip(),
        cell_block.strip(),
        atomic_species_block.strip(),
        formatted_positions.strip(),
        kpoints_block_new.strip()
    ]) + "\n"

    # === Write to output file ===
    basename = os.path.splitext(os.path.basename(input_file))[0]
    out_filename = f"{basename}_SUPER_{sx}_{sy}_{sz}.in"
    with open(out_filename, 'w') as f:
        f.write(output)

    print(f"\n QE input SUCCESSFULLY written to: {out_filename}")
if __name__ == "__main__":
     main()