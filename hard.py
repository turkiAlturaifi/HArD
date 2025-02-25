import argparse
import pandas as pd
import sqlite3
from rdkit import Chem
import sys

DB_FILE = 'hard.db' # database file
# only print these by defualt 
default_cols = ["sigma_het", "homa", "homo", "lumo", "homo_coeff_ipso", "lumo_coeff_ipso",
                "global_electrophilicity", "global_nucleophilicity", "cm5_het",
                "fraction_buried_volume", "sasa_ipso", "sasa", "volume_inside_sas", "length_sterimol_l"]
# with --electronic
electronic_cols = ["sigma_het", "homa", "homo", "lumo", "homo_coeff_ipso", "lumo_coeff_ipso",
                   "chemical_potential", "global_electrophilicity", "global_nucleophilicity",
                   "dispersion_ipso", "dispersion", "npa_ipso", "npa_ipso_carboxylic_acid",
                   "npa_ipso_carboxylate", "npa_het", "npa_het_carboxylic_acid", "npa_het_carboxylate",
                   "hirshfeld_ipso", "hirshfeld_ipso_carboxylic_acid", "hirshfeld_ipso_carboxylate",
                   "hirshfeld_het", "hirshfeld_het_carboxylic_acid", "hirshfeld_het_carboxylate",
                   "cm5_ipso", "cm5_ipso_carboxylic_acid", "cm5_ipso_carboxylate",
                   "cm5_het", "cm5_het_carboxylic_acid", "cm5_het_carboxylate",
                   "total_dipole_moment", "dipole_moment_max", "dipole_moment_median", "dipole_moment_min",
                   "quadrupole_moment_magnitude", "quadrupole_moment_max", "quadrupole_moment_median", "quadrupole_moment_min"]
# --steric
steric_cols = ["distal_volume", "fraction_buried_volume", "sasa_ipso", "sasa", "volume_inside_sas",
               "length_sterimol_l", "min_width_sterimol_b1", "max_width_sterimol_b5",
               "c_ipso-h_bond_length", "c_ipso-c_bond_length_carboxylic_acid", "c_ipso-c_bond_length_carboxylate"]

fingerprint_cols = ["number_of_rings", "number_of_heteroatoms", "number_of_hydrogen_bond_donors",
                    "number_of_hydrogen_bond_acceptors", "topological_polar_surface_area_tpsa", "molecular_weight",
                    "logp", "rotatable_bonds_count", "number_of_heavy_atoms", "number_of_carbon_atoms",
                    "number_of_nitrogen_atoms", "number_of_oxygen_atoms", "number_of_sulfur_atoms",
                    "fraction_of_aromatic_bonds", "fraction_of_sp3_hybridized_carbons", "fraction_of_sp2_hybridized_carbons"]

# --help
parser = argparse.ArgumentParser(
    description="This script prints heteroaryl descriptors from the database using a compound's SMILES string.\n\n"
                "You can search for:\n"
                "  1. A specific regioisomer: To search for a specific regioisomer, use a wildcard in the SMILES string (represented by letter \"A\" in ChemDraw). "
                "For example, for 2-pyridyl, you would use \"[*]C1=NC=CC=C1\".\n\n"
                "  2. All regioisomers: For example, if you input pyridine (\"C1=CC=CN=C1\"), it will return 2-pyridyl, 3-pyridyl, and 4-pyridyl.\n\n"
                "Examples:\n"
                "  python hard.py --smi \"C1=CNC=C1\" # search for all regioisomers\n"
                "  python hard.py --smi \"[*]C1=NC=CC=C1\" --electronic   # 2-pyridyl electronic descriptors\n"
                "  python hard.py --smi \"C1=CC=CN=C1\"   # returns 2-, 3-, and 4-pyridyl\n"
                "  python hard.py --smi \"C1=CNC=C1\" --export output.csv\n",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument('--smi', '-s', type=str, default="C1=CNC=C1", 
                    help="SMILES string to match (default: C1=CNC=C1)")
parser.add_argument('--export', type=str, 
                    help="optional: export full output data to CSV file (provide file name)")
# take smiles as a list from csv or txt file
parser.add_argument('--input_file', '--infile', type=str,
                    help="optional: Provide a CSV or TXT file containing a list of SMILES (one per row). If a header is present, it will be ignored.\n"
                         "a csv file with the descriptors will be automatically generated.")
parser.add_argument('--full', action='store_true', help="full data")
parser.add_argument('--electronic', '--elec', action='store_true', help="only electronic descriptors")
parser.add_argument('--steric', action='store_true', help="only steric descriptors")
parser.add_argument('--fingerprint', action='store_true', help="only fingerprint descriptors")

args = parser.parse_args()

def process_smiles(smiles_str):
    mol = Chem.MolFromSmiles(smiles_str) # convert to canonical smiles
    if mol is None:
        return None 
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    # if the SMILES has a wildcard, then use the 'smiles' column. if not use 'smiles_parent'
    # the latter will give all regioisomers
    wildcard_present = any(atom.GetSymbol() == '*' for atom in mol.GetAtoms())
    column_to_search = "smiles" if wildcard_present else "smiles_parent"
    # SQLite db
    conn = sqlite3.connect(DB_FILE)
    query = f"SELECT * FROM descriptors WHERE {column_to_search} = ?"
    df = pd.read_sql_query(query, conn, params=(canonical_smiles,))
    conn.close()
    if df.empty:
        return None
    return df

def format_df(df, display_cols):
    # only print the requested columns
    display_cols = [col for col in display_cols if col in df.columns]
    df_display = df[display_cols].copy()
    for col in df_display.select_dtypes(include="number").columns:
        if any(kw in col.lower() for kw in ["npa", "cm5", "hirshfeld"]):
            df_display[col] = df_display[col].apply(lambda x: f"{x:.3f}" if pd.notnull(x) else "")
        else:
            df_display[col] = df_display[col].apply(lambda x: f"{x:.2f}" if pd.notnull(x) else "")
    # vertical string for each row
    output_lines = []
    for idx, row in df_display.iterrows():
        output_lines.append(f"Row {idx+1}:")
        for col in df_display.columns:
            output_lines.append(f"  {col}: {row[col]}")
        output_lines.append("-" * 40)
    return "\n".join(output_lines)

# If an input file is provided, process each SMILES from the file
if args.input_file:
    try:
        import os
        # if the file is a .txt file, read it as plain text
        if args.input_file.lower().endswith('.txt'):
            with open(args.input_file, "r") as f:
                smiles_list = [line.strip() for line in f if line.strip()]
        else:
            # otherwise it's csv 
            import csv
            with open(args.input_file, "r", newline="") as f:
                sample = f.read(1024)
                f.seek(0)
                sniffer = csv.Sniffer()
                has_header = sniffer.has_header(sample)
            df_in = pd.read_csv(args.input_file, header=0 if has_header else None)
            # convert the first column to string
            smiles_list = df_in.iloc[:, 0].astype(str).tolist()
    except Exception as e:
        print(f"Error reading input file: {e}")
        exit(1)
    
    # process each SMILES and store results
    results = []
    for smi in smiles_list:
        res_df = process_smiles(smi)
        if res_df is not None:
            # the input smiles, added as a column
            res_df.insert(0, "input_smiles", smi)
            results.append(res_df)
        else:
            # dataframe for just the input SMILES and no descriptors
            results.append(pd.DataFrame({"input_smiles": [smi]}))
    if results:
        final_df = pd.concat(results, ignore_index=True)
    else:
        print("No valid SMILES processed.")
        exit(0)
    
    # full data to csv
    export_filename = args.export if args.export else "output.csv"
    final_df.to_csv(export_filename, index=False)
    print(f"Exported data for {len(smiles_list)} SMILES to {export_filename}")
    exit(0)

# If no input file is provided, process a single SMILES from the command line
# convert to canonical smiles
mol = Chem.MolFromSmiles(args.smi)
if mol is None:
    print("Invalid SMILES string provided.")
    exit(1)
canonical_smiles = Chem.MolToSmiles(mol, canonical=True)

# if the SMILES has a wildcard, then use the 'smiles' column. if not use 'smiles_parent'
# the latter will give all regioisomers
wildcard_present = any(atom.GetSymbol() == '*' for atom in mol.GetAtoms())
column_to_search = "smiles" if wildcard_present else "smiles_parent"

# SQLite db
conn = sqlite3.connect(DB_FILE)
query = f"SELECT * FROM descriptors WHERE {column_to_search} = ?"
matched_df = pd.read_sql_query(query, conn, params=(canonical_smiles,))
conn.close()

if matched_df.empty:
    print("No match found")
    exit(0)

base_cols = ["id", "smiles", "parent"] #always print id, smiles, parent.
if args.full:
    display_cols = matched_df.columns.tolist()  # all columns
elif args.electronic:
    display_cols = base_cols + electronic_cols
elif args.steric:
    display_cols = base_cols + steric_cols
elif args.fingerprint:
    display_cols = base_cols + fingerprint_cols
else:
    display_cols = base_cols + default_cols

display_cols = [col for col in display_cols if col in matched_df.columns]

# format everything as two decmial points unless it's a charge then do three:
df_display = matched_df[display_cols].copy()
for col in df_display.select_dtypes(include="number").columns:
    if any(kw in col.lower() for kw in ["npa", "cm5", "hirshfeld"]):
        df_display[col] = df_display[col].apply(lambda x: f"{x:.3f}" if pd.notnull(x) else "")
    else:
        df_display[col] = df_display[col].apply(lambda x: f"{x:.2f}" if pd.notnull(x) else "")

# print each row vertically in key: value format
for idx, row in df_display.iterrows():
    print(f"Row {idx+1}:")
    for col in df_display.columns:
        print(f"  {col}: {row[col]}")
    print("-" * 40)

# --export to csv
if args.export:
    matched_df.to_csv(args.export, index=False)
    print(f"\nExported full output data to {args.export}")
