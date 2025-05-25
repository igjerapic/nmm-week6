import sys

import numpy as np

def generate_mol(BB_length, SC_length, graft_freq, sequence, filetype='mol'):
    """
    Generate a branched polymer chain and write it to a LAMMPS molecule or data file.

    Parameters:
        BB_length (int): Number of backbone beads.
        SC_length (int): Number of beads in each sidechain.
        graft_freq (int): Frequency of grafting sidechains (e.g., every 5 beads).
        sequence (str): Sequence for backbone types (e.g., 'AB', 'ABA', etc.).
        filetype (str): Either 'mol' or 'data'.
    """

    # --- Constants ---
    l = 0.97  # bond length
    n = BB_length
    s = SC_length
    interval = graft_freq

    # --- Graft sites ---
    graft_sites = np.arange(0, n, interval)
    total_side_beads = s * len(graft_sites) if s > 0 else 0
    total_beads = n + total_side_beads

    # --- Atom positions ---
    pos = np.zeros((total_beads, 3))
    for i in range(n):
        pos[i][0] = l * i

    if s > 0:
        direction = 1
        for i, site in enumerate(graft_sites):
            direction *= -1
            for j in range(s):
                idx = n + i * s + j
                pos[idx][0] = pos[site][0]
                pos[idx][1] = l * (j + 1) * direction

    # --- Atom types ---
    unique_blocks = sorted(set(sequence))
    block_map = {char: idx + 1 for idx, char in enumerate(unique_blocks)}

    seq_len = len(sequence)
    block_size = n // seq_len
    rem = n % seq_len

    # Assign types to backbone
    types = []
    for i in range(seq_len):
        block_type = block_map[sequence[i]]
        count = block_size + (rem if i == seq_len - 1 else 0)
        types += [block_type] * count

    if len(types) != n:
        raise ValueError("Internal error: Backbone type assignment mismatch.")

    # Assign matching sidechain types
    if s > 0:
        for site in graft_sites:
            sc_type = types[site]
            types += [sc_type] * s

    # --- Bond list ---
    bonda = []
    bondb = []

    # Backbone bonds
    for i in range(n - 1):
        bonda.append(i + 1)
        bondb.append(i + 2)

    # Sidechain bonds
    if s > 0:
        for i, site in enumerate(graft_sites):
            graft_idx = site + 1
            start = n + i * s + 1
            bonda.append(graft_idx)
            bondb.append(start)
            for j in range(s - 1):
                bonda.append(start + j)
                bondb.append(start + j + 1)

    # --- Write to file ---
    filename = f'polymer_{sequence}.{filetype}'
    with open(filename, 'w') as f:
        if filetype == 'mol':
            f.write('\n\n')
            f.write(f'{total_beads} atoms\n')
            f.write(f'{len(bonda)} bonds\n\n')
            f.write('Types\n\n')
            for i in range(total_beads):
                f.write(f'{i + 1} {types[i]}\n')

            f.write('\nCoords\n\n')
            for i in range(total_beads):
                f.write(f'{i + 1} {pos[i][0]} {pos[i][1]} {pos[i][2]}\n')

            f.write('\nBonds\n\n')
            for idx, (a, b) in enumerate(zip(bonda, bondb), start=1):
                f.write(f'{idx} 1 {a} {b}\n')

        elif filetype == 'data':
            f.write('LAMMPS Description\n\n')
            f.write(f'{total_beads} atoms\n')
            f.write(f'{len(bonda)} bonds\n')
            f.write(f'{len(set(types))} atom types\n')
            f.write(f'1 bond types\n\n')

            # Box size (simple guess from backbone length)
            xhi = pos[:, 0].max() + 10
            yhi = pos[:, 1].max() + 10
            f.write(f'0.0 {xhi:.2f} xlo xhi\n')
            f.write(f'0.0 {yhi:.2f} ylo yhi\n')
            f.write(f'-5.0 5.0 zlo zhi\n\n')

            f.write('Masses\n\n')
            for t in sorted(set(types)):
                f.write(f'{t} 1.0\n')
            f.write('\n')

            f.write('Atoms # full\n\n')
            for i in range(total_beads):
                f.write(f'{i + 1} 1 {types[i]} 0.0 {pos[i][0]} {pos[i][1]} {pos[i][2]}\n')

            f.write('\nBonds\n\n')
            for idx, (a, b) in enumerate(zip(bonda, bondb), start=1):
                f.write(f'{idx} 1 {a} {b}\n')

        else:
            raise ValueError("filetype must be 'mol' or 'data'")

    print(f'File written: {filename}')


# Example usage
if __name__ == '__main__':
#    generate_mol(BB_length=40, SC_length=5, graft_freq=1, sequence='ABAB', filetype='data')
    chain_length = int(sys.argv[1])
    generate_mol(BB_length=chain_length, SC_length=0, graft_freq=1, sequence='A', filetype='data')
    generate_mol(BB_length=chain_length, SC_length=0, graft_freq=1, sequence='A', filetype='mol')
#    generate_mol(BB_length=40, SC_length=5, graft_freq=1, sequence='ABAB', filetype='mol')
#    generate_mol(BB_length=20, SC_length=0, graft_freq=5, sequence='AB', filetype='data')

