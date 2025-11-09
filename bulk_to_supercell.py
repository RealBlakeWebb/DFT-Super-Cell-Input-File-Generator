"""
Name: Blake Webb
Date: 6/26/2025 -- Updated 11/7/2025
Project: Bulk to Super Cell input generator (Enhanced Multi-ibrav Support)
Description: Program takes bulk input file and outputs scaled supercell input file
             Supports multiple Bravais lattice types (ibrav)
"""

from decimal import Decimal, getcontext
import os
import re

# Precision control for clean decimal formatting
getcontext().prec = 20

# ================================
# === BRAVAIS LATTICE CONFIG ====
# ================================

IBRAV_CONFIG = {
    0: {
        "name": "Free (custom cell)",
        "celldm_params": [],
        "scaling_method": "use_cell_parameters"
    },
    1: {
        "name": "Cubic P (sc)",
        "celldm_params": [1],
        "scaling_method": "simple_cubic"
    },
    2: {
        "name": "Cubic F (fcc)",
        "celldm_params": [1],
        "scaling_method": "cubic"
    },
    3: {
        "name": "Cubic I (bcc)",
        "celldm_params": [1],
        "scaling_method": "cubic"
    },
    4: {
        "name": "Hexagonal/Trigonal P",
        "celldm_params": [1, 3],
        "scaling_method": "hexagonal"
    },
    5: {
        "name": "Trigonal R (3fold axis c)",
        "celldm_params": [1, 4],
        "scaling_method": "trigonal"
    },
    6: {
        "name": "Tetragonal P (st)",
        "celldm_params": [1, 3],
        "scaling_method": "tetragonal"
    },
    7: {
        "name": "Tetragonal I (bct)",
        "celldm_params": [1, 3],
        "scaling_method": "tetragonal"
    },
    8: {
        "name": "Orthorhombic P",
        "celldm_params": [1, 2, 3],
        "scaling_method": "orthorhombic"
    },
    9: {
        "name": "Orthorhombic base-centered",
        "celldm_params": [1, 2, 3],
        "scaling_method": "orthorhombic"
    },
    10: {
        "name": "Orthorhombic face-centered",
        "celldm_params": [1, 2, 3],
        "scaling_method": "orthorhombic"
    },
    11: {
        "name": "Orthorhombic body-centered",
        "celldm_params": [1, 2, 3],
        "scaling_method": "orthorhombic"
    },
    12: {
        "name": "Monoclinic P",
        "celldm_params": [1, 2, 3, 4],
        "scaling_method": "monoclinic"
    },
    13: {
        "name": "Monoclinic base-centered",
        "celldm_params": [1, 2, 3, 4],
        "scaling_method": "monoclinic"
    },
    14: {
        "name": "Triclinic",
        "celldm_params": [1, 2, 3, 4, 5, 6],
        "scaling_method": "triclinic"
    }
}

# ================================
# === HELPER FUNCTIONS ==========
# ================================

def format_fractional_coord(val):
    """Format a fractional or float string to 9 decimal places for QE."""
    if '/' in val:
        num, denom = map(int, val.split('/'))
        value = Decimal(num) / Decimal(denom)
    else:
        value = Decimal(str(val))
    return f"{value:.9f}"

def read_input_file(filepath):
    """Read input file and return contents."""
    with open(filepath, 'r') as f:
        return f.read()

def extract_block(text, block_name):
    """Extract a namelist block from QE input file."""
    pattern = rf"&{block_name}\b(.*?)(?:\n\s*/)"
    match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
    if match:
        return f"&{block_name}{match.group(1)}\n/"
    else:
        print(f"⚠️  Warning: Block &{block_name} not found.")
        return ""

def extract_named_section(text, name):
    """Extract named sections like ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS."""
    pattern = rf"{name}.*?\n((?:.|\n)*?)(?:\n\n|\Z)"
    match = re.search(pattern, text, re.IGNORECASE)
    if match:
        header = re.search(rf"{name}.*?\n", text, re.IGNORECASE).group(0).strip()
        return f"{header}\n{match.group(1).strip()}"
    return ""

def parse_atomic_positions(block):
    """Parse atomic positions from ATOMIC_POSITIONS block."""
    lines = block.strip().splitlines()
    atoms = []
    if not lines or not lines[0].lower().startswith("atomic_positions"):
        print("⚠️  Warning: ATOMIC_POSITIONS block not found or malformed.")
        return []
    for line in lines[1:]:
        if line.strip():
            parts = line.split()
            if len(parts) >= 4:
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            else:
                print(f"⚠️  Skipping malformed atomic line: {line}")
    return atoms

def parse_celldm(block):
    """Parse celldm parameters from SYSTEM block."""
    celldm = {}
    matches = re.findall(r'celldm\((\d+)\)\s*=\s*([0-9.eE+-]+)', block, re.IGNORECASE)
    for idx, val in matches:
        celldm[int(idx)] = float(val)
    return celldm

def parse_ibrav(system_block):
    """Extract ibrav value from SYSTEM block."""
    match = re.search(r'ibrav\s*=\s*(-?\d+)', system_block, re.IGNORECASE)
    if match:
        return int(match.group(1))
    print("⚠️  Warning: ibrav not found in SYSTEM block. Defaulting to 0.")
    return 0

def scale_atomic_positions(atoms, sx, sy, sz):
    """Scale atomic positions for supercell construction."""
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

def format_atomic_positions(atoms):
    """Format atomic positions for output."""
    lines = ["ATOMIC_POSITIONS (crystal)"]
    for label, x, y, z in atoms:
        lines.append(f"{label:<3} {x:.9f} {y:.9f} {z:.9f}")
    return "\n".join(lines)

def update_celldm_hexagonal(celldm1, celldm3, sx, sy, sz):
    """Update celldm for hexagonal/trigonal systems (ibrav=4)."""
    new_celldm1 = celldm1 * sx
    new_celldm3 = celldm3 * sz / sx
    return {1: new_celldm1, 3: new_celldm3}

def update_celldm_cubic(celldm1, sx, sy, sz):
    """Update celldm for cubic systems (ibrav=1,2,3)."""
    # For cubic systems, assume uniform scaling or take maximum
    scale = max(sx, sy, sz) if sx == sy == sz else sx
    new_celldm1 = celldm1 * scale
    return {1: new_celldm1}

def update_celldm_tetragonal(celldm1, celldm3, sx, sy, sz):
    """Update celldm for tetragonal systems (ibrav=6,7)."""
    new_celldm1 = celldm1 * sx
    new_celldm3 = celldm3 * sz / sx
    return {1: new_celldm1, 3: new_celldm3}

def update_celldm_orthorhombic(celldm1, celldm2, celldm3, sx, sy, sz):
    """Update celldm for orthorhombic systems (ibrav=8,9,10,11)."""
    new_celldm1 = celldm1 * sx
    new_celldm2 = celldm2 * sy / sx
    new_celldm3 = celldm3 * sz / sx
    return {1: new_celldm1, 2: new_celldm2, 3: new_celldm3}

def update_celldm_monoclinic(celldm_dict, sx, sy, sz):
    """Update celldm for monoclinic systems (ibrav=12,13)."""
    celldm1 = celldm_dict.get(1, 1.0)
    celldm2 = celldm_dict.get(2, 1.0)
    celldm3 = celldm_dict.get(3, 1.0)
    celldm4 = celldm_dict.get(4, 0.0)  # cos(angle)
    
    new_celldm1 = celldm1 * sx
    new_celldm2 = celldm2 * sy / sx
    new_celldm3 = celldm3 * sz / sx
    
    return {1: new_celldm1, 2: new_celldm2, 3: new_celldm3, 4: celldm4}

def update_celldm_triclinic(celldm_dict, sx, sy, sz):
    """Update celldm for triclinic systems (ibrav=14)."""
    celldm1 = celldm_dict.get(1, 1.0)
    celldm2 = celldm_dict.get(2, 1.0)
    celldm3 = celldm_dict.get(3, 1.0)
    celldm4 = celldm_dict.get(4, 0.0)
    celldm5 = celldm_dict.get(5, 0.0)
    celldm6 = celldm_dict.get(6, 0.0)
    
    new_celldm1 = celldm1 * sx
    new_celldm2 = celldm2 * sy / sx
    new_celldm3 = celldm3 * sz / sx
    
    return {1: new_celldm1, 2: new_celldm2, 3: new_celldm3, 
            4: celldm4, 5: celldm5, 6: celldm6}

def update_celldm_by_ibrav(ibrav, celldm_dict, sx, sy, sz):
    """Route to appropriate celldm update function based on ibrav."""
    if ibrav in [1, 2, 3, -3]:  # Cubic
        celldm1 = celldm_dict.get(1, 1.0)
        return update_celldm_cubic(celldm1, sx, sy, sz)
    
    elif ibrav == 4:  # Hexagonal
        celldm1 = celldm_dict.get(1, 1.0)
        celldm3 = celldm_dict.get(3, 1.0)
        return update_celldm_hexagonal(celldm1, celldm3, sx, sy, sz)
    
    elif ibrav in [6, 7]:  # Tetragonal
        celldm1 = celldm_dict.get(1, 1.0)
        celldm3 = celldm_dict.get(3, 1.0)
        return update_celldm_tetragonal(celldm1, celldm3, sx, sy, sz)
    
    elif ibrav in [8, 9, -9, 91, 10, 11]:  # Orthorhombic
        celldm1 = celldm_dict.get(1, 1.0)
        celldm2 = celldm_dict.get(2, 1.0)
        celldm3 = celldm_dict.get(3, 1.0)
        return update_celldm_orthorhombic(celldm1, celldm2, celldm3, sx, sy, sz)
    
    elif ibrav in [12, -12, 13, -13]:  # Monoclinic
        return update_celldm_monoclinic(celldm_dict, sx, sy, sz)
    
    elif ibrav == 14:  # Triclinic
        return update_celldm_triclinic(celldm_dict, sx, sy, sz)
    
    elif ibrav in [5, -5]:  # Trigonal
        celldm1 = celldm_dict.get(1, 1.0)
        celldm4 = celldm_dict.get(4, 0.0)
        new_celldm1 = celldm1 * max(sx, sy, sz)
        return {1: new_celldm1, 4: celldm4}
    
    else:
        print(f"⚠️  Warning: ibrav={ibrav} not fully supported. Using celldm(1) scaling only.")
        celldm1 = celldm_dict.get(1, 1.0)
        return {1: celldm1 * sx}

def format_celldm_block(celldm_dict):
    """Format celldm parameters as string for output."""
    lines = []
    for idx in sorted(celldm_dict.keys()):
        lines.append(f"  celldm({idx}) = {celldm_dict[idx]},")
    return "\n".join(lines)

def scale_kpoints(kpoints_line, sx, sy, sz):
    """Scale k-points for supercell."""
    k = list(map(int, kpoints_line.split()))
    new_k = [max(1, k[0] // sx), max(1, k[1] // sy), max(1, k[2] // sz)] + k[3:]
    return f"K_POINTS {{automatic}}\n{' '.join(map(str, new_k))}"

def build_system_block(ibrav, celldm_dict, nat, ntyp, system_block_orig):
    """Build new SYSTEM block preserving original parameters."""
    # Extract other parameters from original SYSTEM block
    ecutwfc_match = re.search(r'ecutwfc\s*=\s*([0-9.eE+-]+)', system_block_orig, re.IGNORECASE)
    ecutrho_match = re.search(r'ecutrho\s*=\s*([0-9.eE+-]+)', system_block_orig, re.IGNORECASE)
    occupations_match = re.search(r"occupations\s*=\s*['\"]?(\w+)['\"]?", system_block_orig, re.IGNORECASE)
    smearing_match = re.search(r"smearing\s*=\s*['\"]?(\w+)['\"]?", system_block_orig, re.IGNORECASE)
    degauss_match = re.search(r'degauss\s*=\s*([0-9.eE+-]+)', system_block_orig, re.IGNORECASE)
    input_dft_match = re.search(r"input_dft\s*=\s*['\"]?(\w+)['\"]?", system_block_orig, re.IGNORECASE)
    
    ecutwfc = ecutwfc_match.group(1) if ecutwfc_match else "40"
    ecutrho = ecutrho_match.group(1) if ecutrho_match else "320"
    occupations = occupations_match.group(1) if occupations_match else "smearing"
    smearing = smearing_match.group(1) if smearing_match else "gaussian"
    degauss = degauss_match.group(1) if degauss_match else "0.015"
    input_dft = input_dft_match.group(1) if input_dft_match else "wc"
    
    celldm_str = format_celldm_block(celldm_dict)
    
    system_block = f"""&SYSTEM
  ibrav = {ibrav},
{celldm_str}
  nat = {nat},
  ntyp = {ntyp},
  ecutwfc = {ecutwfc},
  ecutrho = {ecutrho},
  occupations = '{occupations}',
  smearing = '{smearing}',
  degauss = {degauss},
  input_dft = '{input_dft}',
/"""
    
    return system_block

# ================================
# === MAIN WORKFLOW =============
# ================================

def main():
    print("\n" + "="*60)
    print("  Quantum ESPRESSO Supercell Input Generator (Multi-ibrav)")
    print("="*60 + "\n")

    # === Load and read original file ===
    input_file = input("Enter path to QE bulk input file: ").strip()
    if not os.path.exists(input_file):
        print(f"❌ Error: File '{input_file}' not found.")
        return
    
    text = read_input_file(input_file)

    # === Prompt for scaling ===
    scaling = input("Enter scaling [format: X,Y,Z]: ")
    try:
        sx, sy, sz = map(int, scaling.replace(",", " ").split())
        print(f"✓ Scaling factors: {sx}x{sy}x{sz}")
    except:
        print("❌ Error: Invalid scaling format. Use format like '2,2,2' or '2 2 2'")
        return

    # === Extract relevant blocks from original ===
    control_block = extract_block(text, "CONTROL")
    system_block_orig = extract_block(text, "SYSTEM")
    electrons_block = extract_block(text, "ELECTRONS")
    ions_block = extract_block(text, "IONS")
    cell_block = extract_block(text, "CELL")
    atomic_species_block = extract_named_section(text, "ATOMIC_SPECIES")
    atomic_positions_block = extract_named_section(text, "ATOMIC_POSITIONS (crystal)")
    k_points_block = extract_named_section(text, "K_POINTS {automatic}")

    # === Determine ibrav ===
    ibrav = parse_ibrav(system_block_orig)
    if ibrav in IBRAV_CONFIG:
        print(f"✓ Detected ibrav = {ibrav} ({IBRAV_CONFIG[ibrav]['name']})")
    else:
        print(f"⚠️  ibrav = {ibrav} detected but not in standard list")

    # === Parse and scale atomic positions ===
    atoms = parse_atomic_positions(atomic_positions_block)
    if not atoms:
        print("❌ Error: No atoms found. Cannot proceed.")
        return
    
    scaled_atoms = scale_atomic_positions(atoms, sx, sy, sz)
    nat = len(scaled_atoms)
    formatted_positions = format_atomic_positions(scaled_atoms)
    print(f"✓ Scaled {len(atoms)} atoms → {nat} atoms in supercell")

    # === Extract and scale celldm ===
    celldm_dict = parse_celldm(system_block_orig)
    celldm_dict_new = update_celldm_by_ibrav(ibrav, celldm_dict, sx, sy, sz)
    print(f"✓ Updated cell parameters: {celldm_dict_new}")

    # === Get ntyp (number of atomic species) ===
    ntyp = sum(1 for line in atomic_species_block.splitlines() 
               if line.strip() and not line.startswith("ATOMIC_SPECIES"))

    # === Build new SYSTEM block ===
    system_block_new = build_system_block(ibrav, celldm_dict_new, nat, ntyp, system_block_orig)

    # === Adjust and format K_POINTS ===
    if k_points_block:
        kpoints_lines = k_points_block.strip().splitlines()
        if len(kpoints_lines) > 1:
            kpoints_line = kpoints_lines[1]
            kpoints_block_new = scale_kpoints(kpoints_line, sx, sy, sz)
            print(f"✓ K-points scaled")
        else:
            kpoints_block_new = k_points_block
            print("⚠️  K-points block found but cannot scale automatically")
    else:
        kpoints_block_new = ""
        print("⚠️  No K_POINTS block found")

    # === Assemble all blocks ===
    output_parts = [
        control_block.strip(),
        system_block_new.strip(),
        electrons_block.strip(),
        ions_block.strip(),
        cell_block.strip(),
        atomic_species_block.strip(),
        formatted_positions.strip(),
        kpoints_block_new.strip()
    ]
    
    output = "\n\n".join([part for part in output_parts if part]) + "\n"

    # === Write to output file ===
    basename = os.path.splitext(os.path.basename(input_file))[0]
    out_filename = f"{basename}_SUPER_{sx}x{sy}x{sz}.in"
    with open(out_filename, 'w') as f:
        f.write(output)

    print(f"\n✅ QE input SUCCESSFULLY written to: {out_filename}")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
