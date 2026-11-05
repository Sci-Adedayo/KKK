#!/usr/bin/env python3
import os
import re
import glob
import sys
import urllib.request
import pandas as pd
from pymol import cmd

OUTPUT_DIR = "mutants"
PDB_WT = "4zgc.pdb"
os.makedirs(OUTPUT_DIR, exist_ok=True)

AA_MAP = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY',
    'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL',
    'W': 'TRP', 'Y': 'TYR'
}

# -----------------------------------------------------------------------------
# 0. Dynamic CSV File Detection
# -----------------------------------------------------------------------------
csv_files = glob.glob("*.csv")
if not csv_files:
    print("[!] Error: No CSV file found in the current working directory.", file=sys.stderr)
    sys.exit(1)

target_csv = None
mutation_col = None

# Prioritize finding a CSV containing a 'Mutation' column (case-insensitive)
for csv_path in csv_files:
    try:
        df_check = pd.read_csv(csv_path)
        for col in df_check.columns:
            if col.strip().lower() == "mutation":
                target_csv = csv_path
                mutation_col = col
                break
        if target_csv:
            break
    except Exception:
        continue

# Fallback: Use the first available CSV if no explicit 'Mutation' column was matched
if not target_csv:
    target_csv = csv_files[0]
    df_fallback = pd.read_csv(target_csv)
    if not df_fallback.empty:
        mutation_col = df_fallback.columns[0]

print(f"[+] Detected input CSV file: {target_csv} (Target Column: '{mutation_col}')")

# -----------------------------------------------------------------------------
# 1. Download Wild-Type Structure if Missing
# -----------------------------------------------------------------------------
if not os.path.isfile(PDB_WT):
    print(f"[+] Downloading Wild-Type structure ({PDB_WT})...")
    urllib.request.urlretrieve("https://files.rcsb.org/download/4ZGC.pdb", PDB_WT)

# -----------------------------------------------------------------------------
# 2. Prepare and Save Cleaned Isolated Wild-Type (Chain A)
# -----------------------------------------------------------------------------
print("[+] Processing and saving Wild-Type (Chain A)...")
cmd.reinitialize()
cmd.load(PDB_WT, "wt_complex")
cmd.remove("not chain A")
cmd.remove("not polymer.protein")
cmd.remove("not alt '' and not alt 'A'")
cmd.alter("all", "alt=''")
cmd.alter("all", "chain='A'")

wt_out_path = os.path.join(OUTPUT_DIR, "K13_WT_A.pdb")
cmd.save(wt_out_path, "wt_complex")
print(f"[✓] Wild-Type saved to: {wt_out_path}")

# -----------------------------------------------------------------------------
# 3. Re-build Mutants Dynamically from Input CSV
# -----------------------------------------------------------------------------
df = pd.read_csv(target_csv)
if mutation_col not in df.columns:
    print(f"[!] Error: Could not locate a valid mutation column in '{target_csv}'.", file=sys.stderr)
    sys.exit(1)

raw_mutations = df[mutation_col].dropna().astype(str).str.strip().unique().tolist()
print(f"[+] Re-building {len(raw_mutations)} mutants from isolated Chain A...")

successful_mutants = 0
for mut in raw_mutations:
    match = re.match(r"([A-Z])(\d+)([A-Z])", mut)
    if not match:
        continue

    wt_aa = match.group(1)
    res_num = int(match.group(2))
    target_aa = match.group(3)
    target_aa_3 = AA_MAP.get(target_aa)

    if not target_aa_3:
        continue

    cmd.reinitialize()
    cmd.load(PDB_WT, "wt_complex")
    cmd.remove("not chain A")
    cmd.remove("not polymer.protein")
    cmd.remove("not alt '' and not alt 'A'")
    cmd.alter("all", "alt=''")
    cmd.alter("all", "chain='A'")

    cmd.wizard("mutagenesis")
    cmd.refresh_wizard()
    cmd.get_wizard().set_mode(target_aa_3)
    cmd.get_wizard().do_select(f"resi {res_num}")
    cmd.get_wizard().apply()
    cmd.set_wizard()

    out_path = os.path.join(OUTPUT_DIR, f"K13_{wt_aa}{res_num}{target_aa}_A.pdb")
    cmd.save(out_path, "wt_complex")
    successful_mutants += 1

cmd.quit()
print(f"[✓] Completed! Wild-Type and {successful_mutants} mutants generated in '{OUTPUT_DIR}/'.")